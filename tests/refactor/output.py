from numpy import correlate, loadtxt,corrcoef, array_equal


def check_ranked_lists(rank_list_file_python, rank_list_file_r):
    datapy=""
    datar = ""
    fixed_datar = ""
    datapy = loadtxt(rank_list_file_python, dtype = str)
    datar = loadtxt(rank_list_file_r, dtype = str)
    return array_equal(datar, datapy)
    # if (len(datapy) != len(datar) or len(datapy[0]) != len(datar[0])):
    #     return False

    # for i in range(len(datapy)):
    #     if (int(datapy[i][0]) + 1 != int(datar[i][0])): # python indexes are smaller than the indexes of R by 1
    #         return False
    #     if datapy[i][1] != datar[i][1]:
    #         return False
    # return True

def check_correlate(components_file_python, components_file_r):

    datapy = loadtxt(components_file_python, dtype = float)
    datar = loadtxt(components_file_r, dtype = float)
    if (len(datapy) != len(datar) or len(datapy[0]) != len(datar[0])):
         return False
    for i in range(len(datapy[0])):
        cor = corrcoef(datapy[:,i],datar[:,i])
        if (abs( 1- cor[0][1]) > 1e-10):     # or cor[1][0]
            return False #TODO if < 0.95? if == -1?
    
    return True

#python
pc1="/Users/yedidimr/nova/source/glint/k=4t=500refactor.out.components.txt"
pc2="/Users/yedidimr/nova/source/glint/k=5t=500refactor.out.components.txt"
pc3="/Users/yedidimr/nova/source/glint/k=6t=1000refactor.out.components.txt"
pc4="/Users/yedidimr/nova/source/glint/k=6t=500refactor.out.components.txt"
pr1="/Users/yedidimr/nova/source/glint/k=4t=500refactor.out.rankedlist.txt"
pr2="/Users/yedidimr/nova/source/glint/k=5t=500refactor.out.rankedlist.txt"
pr3="/Users/yedidimr/nova/source/glint/k=6t=1000refactor.out.rankedlist.txt"
pr4="/Users/yedidimr/nova/source/glint/k=6t=500refactor.out.rankedlist.txt"
# R
rr1="/Users/yedidimr/nova/source/refactor/refactor.rankedlist.R.k=4.t=500.txt"
rr2="/Users/yedidimr/nova/source/refactor/refactor.rankedlist.R.k=5.t=500.txt"
rr3="/Users/yedidimr/nova/source/refactor/refactor.rankedlist.R.k=6.t=1000.txt"
rr4="/Users/yedidimr/nova/source/refactor/refactor.rankedlist.R.k=6.t=500.txt"
rc1="/Users/yedidimr/nova/source/refactor/refactor.components.R.k=4.t=500.txt"
rc2="/Users/yedidimr/nova/source/refactor/refactor.components.R.k=5.t=500.txt"
rc3="/Users/yedidimr/nova/source/refactor/refactor.components.R.k=6.t=1000.txt"
rc4="/Users/yedidimr/nova/source/refactor/refactor.components.R.k=6.t=500.txt"
# print check_correlate(pc1,rc1)
print check_correlate(pc2,rc2)
# print check_correlate(pc3,rc3)
# print check_correlate(pc4,rc4)
# print check_ranked_lists(pr1,rr1)
print check_ranked_lists(pr2,rr2)
# print check_ranked_lists(pr3,rr3)
# print check_ranked_lists(pr4,rr4)
# print check_correlate("/Users/yedidimr/nova/source/glint/refactor.out.components.txt", "/Users/yedidimr/nova/source/refactor/REUT_R_comp")
# print check_ranked_lists("/Users/yedidimr/nova/source/glint/refactor.out.rankedlist.txt", "/Users/yedidimr/nova/source/refactor/REUT_R_rank_list")