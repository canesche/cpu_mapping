from math import ceil
import traceback
import os


def main(_folder: str, _dest_dir: str):
    files_l = []
    for dir, folder, files in os.walk(_folder):
        for f in files:
            files_l.append([os.path.join(dir, f), f, dir.replace('.','').replace('/','_')])
    # print(files_l)
    stats = []
    for f, n, d in files_l:
        stat = {
            'name': n,
            'bench': f,
            'nodes': 0,
            'edges': 0,
            'ideal_cost': 0,
            'avg_grade': 0,  # output
            'g_hub': {},
            'top_k_hub': [],  # 20%
            'percent_multicast': 0,  # > 2
            'hist_grade': {}  # per grade
        }

        nodes = []
        edges = []
        f = open(f)
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

        # Adding the number of edges and nodes to stats
        stat['nodes'] = len(nodes)
        stat['edges'] = len(edges)
        stat['ideal_cost'] = len(edges)

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
        gr_sum = 0
        hub_sum = 0
        gr_gt = None
        for i in range(len(nodes)):
            nodes[i]['su'] = []
            nodes[i]['pre'] = []
            for j in range(len(edges)):
                if edges[j][1] == nodes[i]['id']:
                    nodes[i]['pre'].append(edges[j][0])
                if edges[j][0] == nodes[i]['id']:
                    nodes[i]['su'].append(edges[j][1])

            if len(nodes[i]['su']) > 2:
                stat['g_hub'][str(i)] = len(nodes[i]['su'])
                hub_sum += 1

            gr_sum += len(nodes[i]['su'])
            if gr_gt is None:
                gr_gt = i
            else:
                if len(nodes[i]['su']) > len(nodes[gr_gt]['su']):
                    gr_gt = i
            # hist:
            if str(len(nodes[i]['su'])) not in stat['hist_grade'].keys():
                stat['hist_grade'][str(len(nodes[i]['su']))] = 1
            else:
                stat['hist_grade'][str(len(nodes[i]['su']))] += 1

        for i in range(int(max(sorted(stat['hist_grade'].keys())))):
            if str(i) not in stat['hist_grade'].keys():
                stat['hist_grade'][str(i)] = 0

        stat['avg_grade'] = gr_sum/stat['nodes']
        stat['percent_multicast'] = hub_sum/stat['nodes']

        top_k_hub = sorted(stat['g_hub'].items(),
                           key=lambda item: item[1], reverse=True)
        k = ceil(len(top_k_hub)*0.2)
        for i in range(k):
            stat['top_k_hub'] = top_k_hub[i]
        stats.append(stat)

    wr_str = ''
    for k in stats[0].keys():
        wr_str += k+';'
    wr_str = wr_str[:-1] + '\n'
    for i in range(len(stats)):
        wr_str += '%s; ' % stats[i]['name']
        wr_str += '%s; ' % stats[i]['bench']
        wr_str += '%d; ' % stats[i]['nodes']
        wr_str += '%d; ' % stats[i]['edges']
        wr_str += '%d; ' % stats[i]['ideal_cost']
        wr_str += '%0.3f; ' % stats[i]['avg_grade']
        wr_str += '%s; ' % str(stats[i]['g_hub'])
        wr_str += '%s; ' % str(stats[i]['top_k_hub'])
        wr_str += '%0.3f; ' % stats[i]['percent_multicast']
        wr_str += '%s \n' % stats[i]['hist_grade']
    wr = open(_dest_dir+'stats_%s.csv'%d, 'w')
    wr.write(wr_str)
    wr.close()

    for j in range(len(stats)):
        wr_str = ''
        for i in range(len(stats[j]['hist_grade'])):
            wr_str += '%d;%d\n' % (i, stats[j]['hist_grade'][str(i)])
        wr = open(_dest_dir+'%s_hist.csv' % stats[j]['name'], 'w')
        wr.write(wr_str)
        wr.close()


if __name__ == '__main__':
    try:
        main('./bench/lisa/dac/', './bench/lisa/stats/')
        a = 1
    except Exception as e:
        print(e)
        traceback.print_exc()
