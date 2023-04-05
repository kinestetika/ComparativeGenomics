import gzip
import re


class TreeNode:
    def __init__(self, id:str='', parent=None, descr='', length = 0):
        self.id = id
        self.descr = descr
        self.parent = parent
        self.length = length
        self.children = []

    def add_to_str(self, my_str, indent):
        my_str += f'{indent}{self.id} - {self.length} - {self.descr}\n'
        for child in self.children:
            my_str = child.add_to_str(my_str, indent+'  ')
        return my_str

    def __str__(self):
        my_str = ''
        my_str = self.add_to_str(my_str, '')
        return my_str


class Tree:
    def __init__(self):
        self.root = TreeNode('root')
        self.node_dict = {'root': self.root}

    def add_node(self, parent, id, descr='') -> TreeNode:
        new_node = TreeNode(id, parent=parent, descr=descr)
        parent.children.append(new_node)
        self.node_dict[id] = new_node
        return new_node

    def __len__(self):
        return len(self.node_dict)


RE_NEWICK_TOKENS = re.compile(r"([():,;'])")
EXPECT_NODE = 'expect_node'
EXPECT_CURRENT_NODE_DESCR = 'expect_current_node_descr'
EXPECT_BRANCH_LENGTH = 'expect_branch_length'

class NewickParser:
    def __init__(self, filename):
        self.filename = filename
        self.handle = None
        self.current_tree = None
        self.current_node = None
        self.status = ''
        self.quoted_str = ''
        self.inside_quotes = False

    def __enter__(self):
        if str(self.filename).endswith('.gz'):
            self.handle = gzip.open(self.filename, 'rt')
        else:
            self.handle = open(self.filename)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.handle.close()

    def __iter__(self):
        for line in self.handle:
            for token in RE_NEWICK_TOKENS.split(line):
                if not token:
                    continue
                #print(f'"{token}"')
                if "'" == token:
                    if self.inside_quotes:
                        if EXPECT_NODE == self.status:
                            self.current_node = self.current_tree.add_node(self.current_node,
                                                                           id=f'node_{len(self.current_tree)}',
                                                                           descr=self.quoted_str)
                        elif EXPECT_CURRENT_NODE_DESCR == self.status:
                            self.current_node.descr = self.quoted_str
                        self.quoted_str = ''
                        self.status = ''
                    else:
                        self.quoted_str = ''
                        self.inside_quotes = False
                elif self.inside_quotes:
                    self.quoted_str += token
                elif ';' == token:
                    final_tree = self.current_tree
                    self.current_tree = None
                    self.current_node = None
                    self.status = ''
                    if final_tree:
                        yield final_tree
                elif '(' == token:
                    if not self.current_tree:
                        self.current_tree = Tree()
                        self.current_node = self.current_tree.root
                    elif self.current_node:
                        self.status = EXPECT_CURRENT_NODE_DESCR
                        self.current_node = self.current_tree.add_node(self.current_node, f'node_{len(self.current_tree)}')
                    self.status = EXPECT_NODE
                elif ')' == token:
                    self.current_node = self.current_node.parent
                    self.status = EXPECT_CURRENT_NODE_DESCR
                elif ':' == token:
                    self.status = EXPECT_BRANCH_LENGTH
                elif ',' == token:
                    self.current_node = self.current_node.parent
                    self.status = EXPECT_NODE
                elif EXPECT_CURRENT_NODE_DESCR == self.status:
                    self.current_node.descr = token
                    self.status = ''
                elif EXPECT_BRANCH_LENGTH == self.status:
                    self.current_node.length = float(token)
                    self.status = ''
                elif EXPECT_NODE == self.status:
                    self.current_node = self.current_tree.add_node(self.current_node, id=f'node_{len(self.current_tree)}',
                                                                   descr=token)
                #if self.current_node:
                #    print(str(self.current_node).strip())
                #else:
                #    print('NO NODE')


def main():
    tree_file = '/home/kinestetika/Desktop/species_tree.nwk'
    with NewickParser(tree_file) as parser:
        for t in parser:
            print(t.root)

if __name__ == "__main__":
    main()