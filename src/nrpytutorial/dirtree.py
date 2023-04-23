import enum
from dataclasses import dataclass
from functools import partial

@enum.unique
class DirtreeNodeKind(enum.Enum):
    FILE = enum.auto()
    DIR = enum.auto()

@dataclass
class DirtreeNode:
    name: str
    kind: DirtreeNodeKind

def make_node(name: str, kind: DirtreeNodeKind) -> DirtreeNode:
    is_kind_valid = kind in DirtreeNodeKind
    if not is_kind_valid:
        raise TypeError('Kind needs be of type NodeKind')
    return DirtreeNode(name, kind)

make_file_node = partial(make_node, kind=DirtreeNodeKind.FILE)
make_dir_node = partial(make_node, kind=DirtreeNodeKind.DIR)
