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
    parent: str

    def to_path(self):
        return NotImplemented

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

    def create_node(self, name, kind, should_backup):
        self._builder.set_name(name)
        self._builder.set_kind(kind)
        self._builder.set_should_backup(should_backup)
        return _builder.build()

    def create_file_node(name, should_backup):
        self._builder.set_name(name)
        self._builder.set_kind(DirtreeNodeKind.FILE)
        self._builder.set_should_backup(should_backup)
        return _builder.build()

    def create_dir_node(name, should_backup):
        self._builder.set_name(name)
        self._builder.set_kind(DirtreeNodeKind.DIR)
        self._builder.set_should_backup(should_backup)
        return _builder.build()

    def create_dir_node(name, should_backup):
        self._builder.set_name(name)
        self._builder.set_kind(DirtreeNodeKind.DIR)
        self._builder.set_should_backup(should_backup)
        return _builder.build()
