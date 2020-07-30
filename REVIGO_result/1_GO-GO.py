"""
提取GO-GO关系
"""
import pandas as pd


def find_related_go(df_data, go_na):
    """
    找名字为go_na的go id和相关GO列表
    :param df_data:treemap数据
    :param go_na:GO名
    :return:GO对应的go id（str），相关GO列表
    """
    inf = df_data[df_data['representative'] == go_na]
    ls = list(inf[inf['description'] == go_na]['term_ID'])
    length = len(ls)
    if length == 0:
        go_id_per_go = 'unknown-'+go_na
    else:
        go_id_per_go = ls[0]
    inf = inf[inf['description'] != go_na]
    return go_id_per_go, list(inf['term_ID'])


def get_similarity(go_id0, go_id1):
    """
    计算两个go间的相似性
    用其对应的基因交集长度比上并集长度
    :param go_id0: go id
    :param go_id1: go id
    :return: 相似性
    """
    set1 = find_go_gene('../info_data/glyma_bp_go_genes.csv', go_id0)
    set2 = find_go_gene('../info_data/glyma_bp_go_genes.csv', go_id1)
    length = len(set1 | set2)
    if length:
        return len(set1 & set2)/len(set1 | set2)
    return 0


def find_go_gene(file_name, go_id_des):
    """
    在f中找GO对应的gene list
    :param file_name: GO与其对应的gene list关系的文件
    :param go_id_des: 要找的GO的id
    :return: go_id对应的gene列表，以集合形式返回
    """
    df_info = pd.DataFrame(pd.read_csv(file_name, sep=','))
    temp = (df_info[df_info['go'] == go_id_des]['gene_ls'])
    if list(temp.values):
        gene_str = list(temp.values)[0]
        set_res = set(eval(gene_str))
    else:
        set_res = set()
    return set_res


def find_go_name(des_id):
    """
    根据go id找到对应的go名字
    :param des_id: go id
    :return: 对应的go名字
    """
    f = 'L7_621_salt_0.001_go_p1e-05_REVIGO_treemap.csv'
    df_info = pd.DataFrame(pd.read_csv(f, sep=','))
    name = (df_info[df_info['term_ID'] == des_id]['description'])
    if len(des_id) == 10:
        return name.values[0]
    # 当go id未知时,des_id为‘unknown-’+go name
    return des_id[8:]


if __name__ == '__main__':
    file = 'L7_621_salt_0.001_go_p1e-05_REVIGO_treemap.csv'
    data = pd.DataFrame(pd.read_csv(file, sep=','))
    # 存放GO-GO关系的表
    GO_interaction_dat = open('salt_GO-GO.csv', 'w')
    GO_interaction_dat.write('GO1,GO2,similarity,type,group_label\n')
    # 存放所有提取的GO节点
    GO_nodes = open('salt_GO_node.csv', 'w')
    GO_nodes.write('node,name,label,group,size(gene_num)\n')
    # core_node_ls  = list(set(list(data['representative'])))
    # 将treemap中的gibberellin biosynthesis、wax biosynthesis替换为
    # gibberellin biosynthetic process和wax biosynthetic process
    core_node_ls = ['cell wall organization', 'gibberellin biosynthetic process',
                    'proteolysis', 'defense response', 'wax biosynthetic process',
                    'iron ion transport', 'cell adhesion', 'positive regulation of development, heterochronic'
                    , 'protein folding']

    # 选出description为核心节点的数据
    dat = data[data['description'].isin(core_node_ls)]
    dat = dat['term_ID,description'.split(',')]
    dat_dic = dat.to_dict()
    # 核心节点的id list
    GO_ids = list(dat['term_ID'])

    for go_name in core_node_ls:
        # 得到每个核心GO的go id和其相关GO
        go_id, related_go_ls = find_related_go(data, go_name)
        # 核心GO写入node文件,去掉related_go_ls为[]的
        if related_go_ls:
            GO_nodes.write('{},{},{},{},{}\n'.format(
                go_id, go_name.replace(',', '/'), go_name.replace(',', '/'), go_id,
                len(find_go_gene('../info_data/glyma_bp_go_genes.csv', go_id))))

        # 有相关GO则写入文件，其余related_go_ls为[]
        for go in related_go_ls:
            # 将上级GO写入文件，有的GO名字含有逗号，用顿号取代
            GO_nodes.write('{},{},{},{},{}\n'.format(
                go, find_go_name(go).replace(',', '/'), 'null',
                go_id, len(find_go_gene('../info_data/glyma_bp_go_genes.csv', go))))

    GO_nodes.close()
    fp = 'salt_GO_node.csv'
    info = pd.DataFrame(pd.read_csv(fp, sep=','))
    # 求node文件中各个GO之间的关系
    for i in range(info.shape[0] - 1):
        for j in range(i + 1, info.shape[0]):
            if info.loc[i, :]['group'] == info.loc[j, :]['group']:
                group = info.loc[j, :]['group']
            else:
                group = 'cross'
            GO_interaction_dat.write('{},{},{},{},{}\n'.format(
                info.loc[i, :]['node'], info.loc[j, :]['node'],
                get_similarity(info.loc[i, :]['node'], info.loc[j, :]['node']),
                'go-go', group))

    GO_interaction_dat.close()
