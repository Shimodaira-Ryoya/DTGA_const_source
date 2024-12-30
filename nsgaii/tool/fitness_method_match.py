"""NSGAII処理において評価値に遺伝子一致率を採用する際に使用します
"""
def list_matcher(a,b):
    """与えられた二つのリストa,bの一致率を計算する

    Args:
        a (list): 比べるリスト
        b (list): 比べるリスト２※二つのリストは同じ長さであること

    Returns:
        match_per: 二つのリストの一致率
    """
    match_point=0#２リストの一致要素数を示す
    for i in range(len(a)):
        if a[i]==b[i]:
            match_point+=1#同要素番号にてa,bの値が一致したら得点
    match_per=match_point/len(a)#一致率の計算
    return  match_per

def match_ind_pop(ind,pop):
    """一つのリストindに対し他の複数のリストpop[other]についてそれぞれ比較し
    indの一致率の平均をもとめる

    Args:
        ind (list): 比較するリスト
        pop (list[list]): 比較する他のリストが入ったリスト(つまり二次元リスト)

    Returns:
        ave: indの一致率平均
    """
    match_per=0
    for other in pop:
        per=list_matcher(ind,other)
        match_per+=per
    match_ave=match_per/len(pop)
    return match_ave