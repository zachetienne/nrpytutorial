import enum
from dataclasses import dataclass
from functools import partial

@enum.unique
class NodeKind(enum.Enum):
    FILE = enum.auto()
    DIR = enum.auto()

@dataclass
class Node:
    name: str
    kind: NodeKind

def make_node(name: str, kind: NodeKind) -> Node:
    is_kind_valid = kind in NodeKind
    if not is_kind_valid:
        raise TypeError('Kind needs be of type NodeKind')
    return Node(name, kind)

make_file_node = partial(make_node, kind=NodeKind.FILE)
