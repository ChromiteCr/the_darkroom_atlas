"""Latent-image formation: the two-site chain and its lookup table (spec D1 2.2)."""

from darkroom_atlas.latent.chain import ChainParams, cluster_distribution, generator
from darkroom_atlas.latent.table import LatentTable, latent_table

__all__ = ["ChainParams", "LatentTable", "cluster_distribution", "generator", "latent_table"]
