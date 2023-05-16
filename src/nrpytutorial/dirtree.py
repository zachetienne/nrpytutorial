from __future__ import annotations

import enum
from dataclasses import dataclass, field, InitVar
from typing import List, Dict
from functools import partial, partialmethod
from pathlib import Path


@enum.unique
class _DirtreeNodeKind(enum.Enum):
    FILE = enum.auto()
    DIR = enum.auto()


@dataclass
class _DirtreeNode:
    name: str
    kind: _DirtreeNodeKind
    should_backup: bool = False

    def to_path(self):
        return Path(self._name)

    def create(self):
        return NotImplemented

def _make_node(name: str,
               kind: _DirtreeNodeKind,
               should_backup: bool = False) -> _DirtreeNode:
    return _DirtreeNode(name, kind, parent, should_backup)

def _make_file_node(name: str, should_backup: bool = False) -> _DirtreeNode:
    return _DirtreeNode(name, _DirtreeNodeKind.FILE, parent, should_backup)

def _make_dir_node(name: str, should_backup: bool = False) -> _DirtreeNode:
    return _DirtreeNode(name, _DirtreeNodeKind.DIR, parent, should_backup)

# @enum.unique
# class DirtreeNodeKind(enum.Enum):
#     FILE = enum.auto()
#     DIR = enum.auto()
#
#
# @enum.unique
# class DirtreeNodeCreateStrategy(enum.Enum):
#     BACKUP = enum.auto()
#     CLEAN = enum.auto()
#     TOUCH = enum.auto()
#
#
# def make_creation_strategy(kind: DirtreeNodeCreateStrategy, node_kind: DirtreeNodeKind):
#     def backup_strategy(path: Path):
#         new_path = path
#         count = 0
#         while new_path.exists():
#             new_path = f"{str(path)}.{count}bak"
#             new_path = Path(new_path)
#             count += 1
#         path.rename(new_path)
#
#     def generic_strategy(path: Path):
#         mthd, config = mthd_lookup[(kind, node_kind)]
#         mthd = getattr(path, mthd)
#         mthd(**config)
#
#     if kind is DirtreeNodeCreateStrategy.BACKUP:
#         return backup_strategy
#
#     mthd_lookup = {
#             (DirtreeNodeCreateStrategy.CLEAN, DirtreeNodeKind.FILE): ("unlink", {"missing_ok": True}),
#             (DirtreeNodeCreateStrategy.CLEAN, DirtreeNodeKind.DIR): ("rmdir", {}),
#             (DirtreeNodeCreateStrategy.TOUCH, DirtreeNodeKind.FILE): ("touch", {"exists_ok": True}),
#             (DirtreeNodeCreateStrategy.TOUCH, DirtreeNodeKind.DIR): ("mkdir", {"exist_ok": True, 'parents': True}),
#             }
#
#     return generic_strategy
#
# @dataclass
# class DirtreeNode:
#     name: str
#     kind: DirtreeNodeKind
#     creation_strategy: DirtreeNodeCreateStrategy
#     _parent: DirtreeNode = field(default_factory=lambda: None)
#     _children: List[DirtreeNode] = field(init=False, default_factory=list)
#
#     def _make_path(self):
#         parent = Path('')
#         parent_head = self._parent
#         while parent_head:
#             parent = Path(parent_head).joinpath(parent)
#             parent_head = parent.parent
#         return parent.joinpath(self.name)
#
#     def create(self):
#         path = self._make_path()
#         do_creation_strategy = make_creation_strategy(self.creation_strategy,
#                                                       self.kind)
#         do_creation_strategy(path)
#
#
#     def update_parent(self, node: DirtreeNode):
#         if node.kind is not DirtreeNodeKind.DIR:
#             raise ValueError('Node parent must be a dir dirtree kind')
#
#         self._parent = node
#
#     @property
#     def parent(self):
#         return self._parent
#
#     def add_child(self, node, *args):
#         if self.kind is not DirtreeNodeKind.DIR:
#             raise ValueError('Node is not dir dirtree kind, it is a leaf')
#         self._children.extend([node, *args])
#
#     @property
#     def children(self):
#         return self._children
#
#
# def make_node(name: str,
#               kind: DirtreeNodeKind,
#               creation_strategy: DirtreeNodeCreateStrategy
#               ) -> DirtreeNode:
#     if kind not in DirtreeNodeKind:
#         raise TypeError(f'Kind {kind} not of type NodeKind')
#
#     if creation_strategy not in DirtreeNodeCreateStrategy:
#         raise TypeError(
#             f'Creation strategy {creation_strategy} not of type NodeKind')
#
#     return DirtreeNode(name, kind, creation_strategy)
#
#
# make_file_node = partial(
#     lambda k, n, s: make_node(
#         n, k, s), DirtreeNodeKind.FILE)
# make_dir_node = partial(
#     lambda k, n, s: make_node(
#         n, k, s), DirtreeNodeKind.DIR)
#
#
# def get_node_factory(kind: DirtreeNodeKind):
#     if kind is DirtreeNodeKind.FILE:
#         return make_file_node
#     return make_dir_node
#
#
# @dataclass
# class Dirtree:
#     _nodes: Dict[str, DirtreeNode] = field(init=False, default_factory=dict)
#     root: InitVar[str]
#
#     def __post_init__(self, root):
#         root = make_dir_node(root, DirtreeNodeCreateStrategy.TOUCH)
#         self.update_root(root)
#
#     def update_root(self, node: DirtreeNode) -> None:
#         if self.test_node_key_used('root'):
#             raise AttributeError(f'root already set')
#
#         if node.kind is not DirtreeNodeKind.DIR:
#             raise TypeError("Root can only be DirtreeNodeKind.DIR")
#
#         if node.parent:
#             raise TypeError('Root should not have parent')
#
#         self._nodes['root'] = node
#
#     def create(self):
#         node_stack = [self._nodes['root']]
#
#         while node_stack:
#             node = node_stack.pop()
#             node.create()
#
#     def create_node(self, key: str, should_backup=False) -> None:
#         if should_backup:
#             self.backup_node(key)
#
#         node = self.get_node(key)
#         path = self.get_node_path(key)
#         if node.kind is DirtreeNodeKind.DIR:
#             path.mkdir(parents=True, exist_ok=True)
#             return
#         path.touch(exist_ok=True)
#         return self.create_node(node)
#
#     def test_node_key_used(self, key: str) -> bool:
#         return key in self._nodes
#
#     def add_node(self, node: DirtreeNode, parent: str, key: str = None):
#         if 'root' not in self._nodes:
#             raise AttributeError('Missing root node, set it first')
#
#         key = key if key else node.name
#         if self.test_node_key_used(key):
#             raise AttributeError(f'Node key {key} already used')
#
#         if not self.test_node_key_used(parent):
#             raise ValueError(f'Parent {parent} node not found')
#
#         parent = self.get_node(parent)
#         node.update_parent(parent)
#         self._nodes[key] = node
#
#     def get_node(self, key: str) -> DirtreeNode:
#         if not self.test_node_key_used(key):
#             raise AttributeError(f'Node key {key} not found')
#         return self._nodes[key]
#
#     get_root = partialmethod(get_node, 'root')
#
#     def get_node_path(self, key: str) -> Path:
#         node = self.get_node(key)
#         node_path = Path(node.name)
#         while node.parent is not None:
#             node_path = Path(node.parent.name).joinpath(node_path)
#             node = node.parent
#         return node_path
#
#     def get_node_str(self, key: str) -> str:
#         return str(self.get_node_path(key))
#
#     def backup_node(self, key: str) -> None:
#         path = self.get_node_path(key)
#         new_path = path
#         count = 0
#         while new_path.exists():
#             new_path = f"{str(path)}.{count}bak"
#             new_path = Path(new_path)
#             count += 1
#         path.rename(new_path)
#         return
#
#     create_root = partialmethod(create_node, 'root')
#
#     def build_node(self,
#                    kind: DirtreeNodeKind,
#                    name: str,
#                    parent: str,
#                    key: str):
#         node_factory = get_node_factory(kind)
#         node = node_factory(name)
#         self.add_node(node, parent, key)
#         self.create_node(key)
