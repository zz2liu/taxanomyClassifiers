"""p0_compiler.py: compile p0 language to c.
"""
class Node(object):
    pass

class Module(Node):
    def __init__(self, doc, node):
        self.doc = doc
        self.node = node


class Stmt(Node):
    def __init__(self, nodes):
        pass

class Printnl(Node):
    pass


def num_nodes(ast):
    if isinstance(ast, Module):
        return 1 + num_nodes(ast.node)
    elif isinstance(ast, Stmt):
        return 1 + sum(num_nodes(c) for c in ast.nodes)
