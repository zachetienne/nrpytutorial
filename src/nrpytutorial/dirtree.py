from __future__ import annotations

import enum
from dataclasses import dataclass, field, InitVar
from typing import List, Dict
from pathlib import Path


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

    def to_path(self):
        parent_node = self.parent
        parent = Path('')
        while parent_node:
            parent = Path(parent_node.data.name).joinpath(parent)
            parent_node = parent_node.parent
        return parent.joinpath(self.data.name)

    def to_string(self):
        path = self.to_path()
        return str(path)

    def create(self):
        return NotImplemented

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
        self._builder.set_name(name)
        self._builder.set_kind(DirtreeNodeKind.FILE)
        self._builder.set_parent(parent)
        self._builder.set_should_backup(should_backup)
        return self._builder.build()

    def construct_dir_node(self,
                        name: str,
                        parent: DirtreeNode = None,
                        should_backup: bool = False) -> DirtreeNode:
        self._builder.set_name(name)
        self._builder.set_kind(DirtreeNodeKind.DIR)
        self._builder.set_parent(parent)
        self._builder.set_should_backup(should_backup)
        return self._builder.build()
