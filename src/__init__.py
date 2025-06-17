"""
AI Knowledge Graph Module.

This module serves as the entry point for the AI Knowledge Graph application,
coordinating data extraction, processing, and storage.
"""
__all__ = ['data_pipeline', 'knowledge_graph_builder']

from . import data_pipeline
from . import knowledge_graph_builder

__version__ = "0.1.0"
