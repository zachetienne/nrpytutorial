from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import List, Dict
from functools import partial, partialmethod
from pathlib import Path


@enum.unique
class DirtreeNodeKind(enum.Enum):
    FILE = enum.auto()
    DIR = enum.auto()

@dataclass
class DirtreeNode:
    name: str
    kind: DirtreeNodeKind
    _parent: DirtreeNode = field(default_factory=lambda: None)
    _children: List[DirtreeNode] = field(init=False, default_factory=list)

    def update_parent(self, node: DirtreeNode):
        if node.kind is not DirtreeNodeKind.DIR:
            raise ValueError('Node parent must be a dir dirtree kind')

        self._parent = node

    @property
    def parent(self):
        return self._parent

    def add_child(self, node, *args):
        if self.kind is not DirtreeNodeKind.DIR:
            raise ValueError('Node is not dir dirtree kind, it is a leaf')
        self._children.extend([node, *args])

    @property
    def children(self):
        return self._children

def make_node(name: str,
              kind: DirtreeNodeKind,
              ) -> DirtreeNode:
    is_kind_valid = kind in DirtreeNodeKind
    if not is_kind_valid:
        raise TypeError('Kind needs be of type NodeKind')
    return DirtreeNode(name, kind)

make_file_node = partial(make_node, kind=DirtreeNodeKind.FILE)
make_dir_node = partial(make_node, kind=DirtreeNodeKind.DIR)

@dataclass
class Dirtree:
    _nodes: Dict[str, DirtreeNode] = field(init=False, default_factory=dict)

    def test_node_key_used(self, key: str) -> bool:
        return key in self._nodes

    def set_root(self, node: DirtreeNode) -> None:
        if self.test_node_key_used('root'):
            raise AttributeError(f'root already set')

        if node.kind is not DirtreeNodeKind.DIR:
            raise TypeError("Root can only be DirtreeNodeKind.DIR")

        if node.parent:
            raise TypeError('Root should not have parent')

        self._nodes['root'] = node

    def add_node(self, node: DirtreeNode, parent: str, key: str = None):
        if not 'root' in self._nodes:
            raise AttributeError('Missing root node, set it first')

        key = key if key else node.name
        if self.test_node_key_used(key):
            raise AttributeError(f'Node key {key} already used')

        if not self.test_node_key_used(parent):
            raise ValueError(f'Parent {parent} node not found')

        parent = self.get_node(parent)
        node.update_parent(parent)
        self._nodes[key] = node

    def get_node(self, key: str) -> DirtreeNode:
        if not self.test_node_key_used(key):
            raise AttributeError(f'Node key {key} not found')
        return self._nodes[key]

    get_root = partialmethod(get_node, 'root')

    def get_node_path(self, key: str) -> Path:
        node = self.get_node(key)
        node_path = Path(node.name)
        while node.parent is not None:
            node_path = Path(node.parent.name).joinpath(node_path)
            node = node.parent
        return node_path
