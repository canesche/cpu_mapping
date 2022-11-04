import traceback


def main(dot: str):
    name = 'dot'
    nodes = []
    edges = []
    f = open(dot)
    lines = f.readlines()
    name = lines[0].split(' ')[1]
    lines.pop(0)
    lines.pop(-1)
    nodes_counter = 0
    for i in range(len(lines)):
        lines[i] = lines[i].replace('\n', '').replace(
            ' ', '').replace(';', '').lower()
        p = lines[i].find('[')
        lines[i] = lines[i][:p]

    n_aux = []
    for l in lines:
        if '//' not in l and l != '':
            if '->' not in l:
                nodes.append({'id': str(nodes_counter), 'name': l})
                n_aux.append(l)
                nodes_counter += 1
            else:
                edges.append(l)

    # cleanning nodes and edges
    # edges
    for i in range(len(edges)):
        # spliting the edes in a vector with the name between spaces
        edges[i] = edges[i].split('->')
        # adding orphan nodes to nodes' dictionary

        if edges[i][0] not in n_aux:
            nodes.append({'id': str(nodes_counter),
                         'name': edges[i][0].replace(' ', '')})
            n_aux.append(edges[i][0])
            nodes_counter += 1
        if edges[i][1] not in n_aux:
            nodes.append({'id': str(nodes_counter),
                         'name': edges[i][1].replace(' ', '')})
            n_aux.append(edges[i][1])
            nodes_counter += 1

        #edges[i][0] = ' ' + edges[i][0] + ' '
        #edges[i][1] = ' ' + edges[i][1] + ' '

        ch = 0
        for n in nodes:
            if edges[i][0] == n['name']:
                edges[i][0] = n['id']
                ch += 1
            elif edges[i][1] == n['name']:
                edges[i][1] = n['id']
                ch += 1
            if ch == 2:
                break

    # finding the successors
    for i in range(len(nodes)):
        nodes[i]['su'] = []
        for j in range(len(edges)):
            if edges[j][1] == nodes[i]['id']:
                nodes[i]['su'].append(edges[j][0])

    for i in range(len(nodes)):
        # Finding the type of the node
        if 'mul' in nodes[i]['name']:
            nodes[i]['op'] = 'mul'
        elif 'add' in nodes[i]['name']:
            nodes[i]['op'] = 'add'
        elif 'sub' in nodes[i]['name']:
            nodes[i]['op'] = 'sub'
        elif 'reg' in nodes[i]['name'] or 'lod' in nodes[i]['name'] or 'load' in nodes[i]['name']:
            nodes[i]['op'] = 'reg'
            nodes[i]['name'] = nodes[i]['name'].replace('lod', 'reg')
            nodes[i]['name'] = nodes[i]['name'].replace('load', 'reg')
        elif 'const' in nodes[i]['name'] or 'in' in nodes[i]['name']:
            nodes[i]['op'] = 'in'
            nodes[i]['name'] = nodes[i]['name'].replace('const', 'in')
        elif 'out' in nodes[i]['name'] or 'str' in nodes[i]['name'] or 'store' in nodes[i]['name']:
            nodes[i]['op'] = 'out'
            nodes[i]['name'] = nodes[i]['name'].replace('str', 'out')
            nodes[i]['name'] = nodes[i]['name'].replace('store', 'out')
        else:
            nodes[i]['op'] = 'NONE'

    # finding which nodes are imadiate ondes
    for i in range(len(nodes)):
        nodes[i]['im'] = False
        if nodes[i]['op'] in ['mul', 'add', 'sub'] and len(nodes[i]['su']) < 2:
            nodes[i]['op'] = nodes[i]['op']+'i'
            nodes[i]['im'] = True
    # create the dot string
    dot_str = 'digraph %s {\n' % name
    # first the nodes
    for n in nodes:
        if n['im']:
            dot_str += '%s [label=%s op=%s value=2];\n' % (
                n['id'], n['name'], n['op'])
        else:
            dot_str += '%s [label=%s op=%s];\n' % (
                n['id'], n['name'], n['op'])
    dot_str += '\n'
    # now the edges

    for n in nodes:
        flag = False
        for s in n['su']:
            if flag:
                dot_str += '%s -> %s [port=1 weight=0];\n' % (s, n['id'])
            else:
                dot_str += '%s -> %s [port=0 weight=0];\n' % (s, n['id'])
                flag = True
    dot_str += '}'
    print(dot_str)


if __name__ == '__main__':
    try:
        main('./bench/m_bench/h2v2_smooth.dot')
    except Exception as e:
        print(e)
        traceback.print_exc()
