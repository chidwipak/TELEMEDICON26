"""D-MPNN encoder and Siamese relative-affinity network (extracted unchanged from the internship code, models.py)."""
from __future__ import annotations

import torch
import torch.nn as nn
from torch_geometric.utils import scatter
from torch_geometric.nn import global_mean_pool


class DMPNNEncoder(nn.Module):
    """Directed MPNN. Returns a graph-level embedding [B, hidden]."""

    def __init__(self, node_dim: int, edge_dim: int, hidden: int = 128,
                 depth: int = 4, dropout: float = 0.1):
        super().__init__()
        self.depth = depth
        self.W_i = nn.Linear(node_dim + edge_dim, hidden)
        self.W_h = nn.Linear(hidden, hidden)
        self.W_o = nn.Linear(node_dim + hidden, hidden)
        self.dropout = nn.Dropout(dropout)
        self.act = nn.ReLU()

    @staticmethod
    def _reverse_index(E: int, device) -> torch.Tensor:
        # edges are built as consecutive (i->j, j->i) pairs => reverse = idx XOR 1
        idx = torch.arange(E, device=device)
        return idx ^ 1

    def forward(self, x, edge_index, edge_attr, batch, return_nodes=False):
        src, dst = edge_index[0], edge_index[1]
        E = edge_attr.size(0)
        N = x.size(0)
        rev = self._reverse_index(E, x.device)

        # initial directed-edge hidden state
        h0 = self.act(self.W_i(torch.cat([x[src], edge_attr], dim=1)))  # [E, H]
        h = h0
        for _ in range(self.depth - 1):
            agg = scatter(h, dst, dim=0, dim_size=N, reduce="sum")       # [N, H]
            m = agg[src] - h[rev]                                        # exclude reverse
            h = self.act(h0 + self.W_h(m))
            h = self.dropout(h)

        agg = scatter(h, dst, dim=0, dim_size=N, reduce="sum")           # [N, H]
        h_node = self.act(self.W_o(torch.cat([x, agg], dim=1)))          # [N, H]
        h_node = self.dropout(h_node)
        if return_nodes:
            return h_node                                               # [N, H]
        return global_mean_pool(h_node, batch)


class SiameseAVES(nn.Module):
    """Siamese fragment network predicting relative ddG between two ligands."""

    def __init__(self, node_dim, edge_dim, hidden=128, depth=4, dropout=0.1):
        super().__init__()
        self.encoder = DMPNNEncoder(node_dim, edge_dim, hidden, depth, dropout)
        self.head = nn.Sequential(
            nn.Linear(hidden * 3, hidden), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(hidden, hidden // 2), nn.ReLU(),
            nn.Linear(hidden // 2, 1),
        )

    def encode(self, data):
        return self.encoder(data.x, data.edge_index, data.edge_attr, data.batch)

    def forward(self, data_cand, data_ref):
        g_c = self.encode(data_cand)
        g_r = self.encode(data_ref)
        combined = torch.cat([g_c, g_r, torch.abs(g_c - g_r)], dim=1)
        return self.head(combined).squeeze(-1)
