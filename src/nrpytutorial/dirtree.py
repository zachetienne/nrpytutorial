from __future__ import annotations
import enum
from dataclasses import dataclass, field
from typing import List, Dict
from functools import partial
from pathlib import Path


@enum.unique
class DirtreeNodeKind(enum.Enum):
    FILE = enum.auto()
    DIR = enum.auto()

@dataclass
class DirtreeNode:
    name: str
    kind: DirtreeNodeKind
    parent: DirtreeNode = field(init=False, default_factory=lambda: None)
    children: List[DirtreeNode] = field(init=False, default_factory=list)

def make_node(name: str, kind: DirtreeNodeKind) -> DirtreeNode:
    is_kind_valid = kind in DirtreeNodeKind
    if not is_kind_valid:
        raise TypeError('Kind needs be of type NodeKind')
    return DirtreeNode(name, kind)

make_file_node = partial(make_node, kind=DirtreeNodeKind.FILE)
make_dir_node = partial(make_node, kind=DirtreeNodeKind.DIR)

@dataclass
class Dirtree:
    _nodes: Dict[str, DirtreeNode] = field(init=False, default_factory=dict)

    def _test_node_key_uniqness(self, alias: str) -> bool:
        return alias in self._nodes

    def assure_alias_unique(self, alias: str) -> None:
        is_alias_used = self._test_node_key_uniqness(alias)
        if is_alias_used:
            raise ValueError(f'Node alias {alias} already used')

    def assure_parent_exists(self, parent: str) -> None:
        is_parent_existing = self._test_node_key_uniqness(parent)
        if not is_parent_existing:
            raise ValueError(f'Parent {parent} not found')

    def assure_parent_kind_dir(self, parent: str) -> None:
        if self._nodes[parent].kind is not DirtreeNodeKind.DIR:
            raise TypeError(f'Parent {parent} is not directory')

    def set_root(self, node: DirtreeNode) -> None:
        if 'root' in self._nodes:
            raise AttributeError(f'root already set')

        if node.kind is not DirtreeNodeKind.DIR:
            raise TypeError(f"Root can only be directory")

        self._nodes['root'] = node

    def add_node(self, node: DirtreeNode, parent: str, alias: str = None):
        if not 'root' in self._nodes:
            raise AttributeError('Missing root node, set it first')

        alias = alias if alias else node.name
        self.assure_alias_unique(alias)
        self.assure_parent_exists(parent)
        self.assure_parent_kind_dir(parent)

    def get_node(self, alias: str) -> Path:
        if not alias in self._nodes:
            raise AttributeError(f'No node alias found with alias = {alias}')
        node = self._nodes[alias]
        return node
