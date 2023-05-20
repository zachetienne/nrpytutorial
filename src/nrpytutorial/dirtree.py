from __future__ import annotations

import enum
from dataclasses import dataclass, field, InitVar
from typing import List, Dict
from pathlib import Path
from functools import partialmethod

@enum.unique
class DirtreeNodeKind(enum.Enum):
    FILE = enum.auto()
    DIR = enum.auto()

@dataclass
class DirtreeNodeData:
    name: str
    kind: DirtreeNodeKind
    should_backup: bool = False

@dataclass
class DirtreeNode:
    data: DirtreeNodeData
    parent: DirtreeNode

    @property
    def path(self):
        parent_node = self.parent
        parent = Path('')
        while parent_node:
            parent = Path(parent_node.data.name).joinpath(parent)
            parent_node = parent_node.parent
        return parent.joinpath(self.data.name)

    def __str__(self):
        return str(self.path)

    def backup(self, is_forced=False):
        should_backup = is_forced if is_forced else self.data.should_backup

        if not should_backup:
            return

        count = 0
        backup_template = "{}_bak{}"
        backup_name = backup_template.format(self.path, count)
        while Path(backup_name).exists():
            count += 1
            backup_name = backup_template.format(self.path, count)

        self.path.rename(backup_name)
        return

    def clean(self):
        if not self.path.exists():
            return

        if self.data.kind is DirtreeNodeKind.DIR:
            return

        self.path.unlink()

    def create(self):
        if self.path.exists():
            self.backup()
            self.clean()

        if self.data.kind is DirtreeNodeKind.DIR:
            self.path.mkdir(parents=True, exist_ok=True)
            return

        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.touch()
        return

    @staticmethod
    def get_builder():
        return DirtreeNodeBuilder()

@dataclass
class DirtreeNodeBuilder:

    def set_name(self, name):
        self._name = name

    def set_kind(self, kind):
        self._kind = kind

    def set_should_backup(self, should_backup):
        self._should_backup = should_backup

    def set_parent(self, parent):
        self._parent = parent

    def build(self):
        data = DirtreeNodeData(self._name, self._kind, self._should_backup)
        return DirtreeNode(data, self._parent)

@dataclass
class DirtreeNodeDirector:
    _builder: DirtreeNodeBuilder = field(default_factory=DirtreeNode.get_builder)

    def construct_node(self,
                    name: str,
                    kind: DirtreeNodeKind,
                    parent: DirtreeNode = None,
                    should_backup: bool =False) -> DirtreeNode:
        self._builder.set_name(name)
        self._builder.set_kind(kind)
        self._builder.set_parent(parent)
        self._builder.set_should_backup(should_backup)
        return self._builder.build()

    def construct_file_node(self,
                         name: str,
                         parent: DirtreeNode = None,
                         should_backup: bool = False) -> DirtreeNode:
        node = self.construct_node(name, DirtreeNodeKind.FILE, parent, should_backup)
        return node

    def construct_dir_node(self,
                         name: str,
                         parent: DirtreeNode = None,
                         should_backup: bool = False) -> DirtreeNode:
        node = self.construct_node(name, DirtreeNodeKind.DIR, parent, should_backup)
        return node

@dataclass
class Dirtree:
    _root: DirtreeNode
    _nodes: dict[str, DirtreeNode] = field(init=False, default_factory=dict)

    def __post_init__(self):
        self._nodes['root'] = self._root
    
    def get(self, key: str):
        if not key in self._nodes:
            raise KeyError(f'Node with key: {key} not found')
        return self._nodes.get(key)

    def add(self, node: DirtreeNode, parent_key: str, key: str = None):
        key = key if key else node.data.name
        node.parent = self.get(parent_key)
        self._nodes[key] = node

    add_codes = partialmethod(add, parent_key='root', key='codes')
    add_build = partialmethod(add, parent_key='root', key='build')
    add_bin = partialmethod(add, parent_key='root', key='bin')
    add_output = partialmethod(add, parent_key='root', key='output')

    def create(self):
        node_keys = sorted(self._nodes, key=lambda k: len(self._nodes[k].path.parents), reverse=True)
        nodes = [ self._nodes[k] for k in node_keys ]
        for node in nodes:
            if node.data.kind is not DirtreeNodeKind.DIR:
                continue
            node.create()

        for node in nodes:
            if node.data.kind is not DirtreeNodeKind.FILE:
                continue
            node.create()
