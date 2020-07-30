"""
提取GO和miRNA对应关系
"""
import math
import pandas as pd


def get_from_go_info(info, cal, threshold):
    """
    对于单一模块，从go_info中卡阈值，与treemap中GO取交集，过滤出需要的GO list
    :param info: 从整个文件提取的数据
    :param cal: 计数，第几个模块
    :param threshold: GO的bh_pVal应满足的阈值
    :return: 过滤结果，GO list和每个GO对应的bh_pVal
    """
    go_info = info['go_info']
    # 先取出一个模块，一共27个模块
    go_info_dict = eval(go_info[cal])
    go = []
    # 将 GO ID：bh_pVal 存入字典
    go_bh_dict = {}
    bh = []
    # 过滤出bh_pVal小于阈值的GO
    for val in go_info_dict.values():
        if val['bh_pVal'] < threshold:
            go.append(val['GO_id'])
            go_bh_dict[val['GO_id']] = val['bh_pVal']
    fp2 = '../REVIGO_result/L7_621_salt_0.001_go_p1e-05_REVIGO_treemap.csv'
    data2 = pd.DataFrame(pd.read_csv(fp2, sep=','))
    go_term_id = set(data2['term_ID'])
    # 取得go_info和treemap中GO的交集
    go_ls = go_term_id & set(go)
    for per in go_ls:
        # 从字典中取出GO对应的bh_pVal
        bh.append(go_bh_dict[per])
    return list(go_ls), bh


def find_go_gene(df_info, go_id):
    """
    在df_info中找GO对应的gene list
    :param df_info: GO与其对应的gene list关系的文件数据
    :param go_id: 要找的GO的id
    :return: go_id对应的gene列表，以集合形式返回
    """
    temp = (df_info[df_info['go'] == go_id]['gene_ls'])
    if list(temp.values):
        gene_str = list(temp.values)[0]
        set_res = set(eval(gene_str))
    else:
        set_res = set()
    return set_res


def find_gene_mirna(info, ge_ls):
    """
    从info中找到ge_ls中所有gene对应的所有miRNA
    :param info: gene、miRNA靶向文件
    :param ge_ls: GO列表
    :return: 对应的miRNA集合
    """
    mirna_set = set()
    for ge in ge_ls:
        temp = info[info['Ensembl_ID'] == ge]['miRNA_id']
        mirna_set.update(list(temp))
    return mirna_set


if __name__ == '__main__':
    fp = 'L7_621_salt_stress_p_add_miR_ratio_salt_stress_fdr_p_mark_final_0.001.csv'
    data = pd.DataFrame(pd.read_csv(fp, sep=','))
    fp0 = '../info_data/glyma_bp_go_genes.csv'
    data0 = pd.DataFrame(pd.read_csv(fp0, sep=','))
    fp1 = '../info_data/AA gma2 data labeled final_add_name.csv'
    data_gene_mirna = pd.DataFrame(pd.read_csv(fp1, sep=','))
    fina_file = 'salt_gene_miRNA.csv'
    gene_mirna = open(fina_file, 'w')
    gene_mirna.write('mirna,go,similarity\n')
    for mod in range(27):
        # 每个模块对应的gene集合
        gene_ls_per_mod = set(eval(data['gene_ls'][mod]))
        # 每个模块对应的mirna集合
        mirna_ls_per_mod = set(eval(data['miRNA_ls'][mod]))

        # 获得一个模块中所有过滤后的go list和对应的bh_pVal
        go_list_per_mod, bh_pval = get_from_go_info(data, mod, 0.001)

        # 所有模块中的每个GO对应的gene列表
        mirna_ls = []
        for i, per_go in enumerate(go_list_per_mod):
            # 在glyma_bp_go_genes.csv中得到的，每个GO对应的gene集合
            gene_set_per_go = find_go_gene(data0, per_go)
            # 一个模块中的gene list和每个GO对应的gene list取交集
            gene = list(gene_ls_per_mod & gene_set_per_go)
            # 每个GO对应的所有miRNA集合
            mirna_per_go = find_gene_mirna(data_gene_mirna, gene)
            mirna = list(mirna_ls_per_mod & mirna_per_go)
            mirna_ls.append(mirna)
        for i in range(len(go_list_per_mod)):
            for j in range(len(mirna_ls[i])):
                gene_mirna.write('{},{},{}\n'.format
                                 (mirna_ls[i][j], go_list_per_mod[i], -math.log(bh_pval[i])))
    # 及时关闭，否则会出错
    gene_mirna.close()
    df = pd.DataFrame(pd.read_csv(fina_file, sep=','))
    # 过滤所有列，p值不同，无重复，668行；过滤前两列，剩余614行
    df.drop_duplicates(subset=['mirna', 'go'], inplace=True)

    gene_mirna_new = open('salt_gene_miRNA_new1.csv', 'w')
    gene_mirna_new.write('mirna,go,similarity\n')
    fina_file_new = 'salt_gene_miRNA_new1.csv'
    df.to_csv(fina_file_new, header=True, index=False)
    print(df.shape)

    df.drop_duplicates(subset=['go'], inplace=True)
    print(df.shape)

    gene_mirna_new.close()
