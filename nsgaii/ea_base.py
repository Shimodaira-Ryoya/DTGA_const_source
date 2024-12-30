import random, sys
import copy
import pandas as pd

class Individual(list):
    """一つの遺伝子"""
    def __init__(self, nobj, nvar, vlow, vhigh, initype):
        self.nvar = nvar#評価関数の個数
        self.nobj = nobj#遺伝子の変数の数
        self.fitness = ()#評価関数の値
        self.values=[]#outputするパラメータ
        self.dtinfo=None#個体のdtinfoを保持
        self.rank = [-1.0] * 2#生存選択における個体のランクを示す
        self.penalty=0 #ペナルティ0より大きい場合、実行不可能解に分類する
        if initype == 'binary':#遺伝子の変数タイプ binaryなら変数の範囲は0と1のみ
            for i in range(nvar):
                if random.random() < vhigh :#vhighは初期遺伝子生成においてある変数で1になる確率
                    self.append(1)
                else:
                    self.append(0)
        else:
            sys.exit("Init Error: init " + str(initype))

    def printind(self):
        output = ' '.join(map(str,self.values))#パラメータ
        output = output + ' ' + ' '.join(map(str,self.rank))#非支配ランク，混雑距離
        output = output + ' gene ' + ' '.join(map(str,self))#遺伝子
        return output

class Population(list):
    """遺伝子集団（individualの集団）"""
    def __init__(self, size=0, nobj=0, nvar=0, vlow=0, vhigh=0, initype=0):
        """遺伝子individualをsize分生成"""
        for i in range(size):#sizeは遺伝子の数を示す。
            self.append(Individual(nobj, nvar, vlow, vhigh, initype))
    
    def printpop(self):
        """遺伝子それぞれのパラメータを出力"""
        for i in range(len(self)):
            print(i, self[i].printind())
        
    def fprintbest(self,f, i):
        """遺伝子集団のうちi番目の遺伝子をファイルに書き込み
        Args:
            f : openしたファイルのポインタ
            i : 遺伝子の番号
        """
        f.write(self[i].printind()+'\n')
    
    def fprintpop(self,f):
        """遺伝子集団全部書き込み
        Args:
            f : openしたファイルのポインタ
        """
        for i in range(len(self)):
            f.write(self[i].printind()+'\n')

    def fprintpopcsv(self,pas,columns):
        """遺伝子集団全部書き込み(csvファイル)
        Args:
            pas: 保存するファイル名＋保存するディレクトリ
            columns:valuesの項目名リスト
        """
        #集団個々のパラメータを一つのリストにまとめ、それらを要素として二次元リストdataを作成
        data=[]
        for i in range(len(self)):
            data.append(self[i].values+self[i].rank+["gene"]+self[i])#values+rank+gene
        #遺伝子ビット一つ一つのカラム名を取得
        genename=[]
        for j in range(len(self[0])):
            genename.append("gene["+str(j)+"]")
        #dataをdf化、カラム名を定義fitnessの返り値columnsに以下要素を加えたものとする
        df = pd.DataFrame(data, columns=columns+["f4(match)","front","CD","gene"]+genename)
        #csv化
        df.to_csv(pas, index=False, encoding="utf-8")
        


"""親選択、交叉、突然変異のメソッド"""

"""親選択"""
def binary_tournament(pop, n):
    """遺伝子集団のうち個体をランダムに二つ選び、性能を比較し(トーナメント)、優秀な方を親として返す
        fitnessで比較を行う
    Args:
        pop (gene.population): 遺伝子集団
        n (int): トーナメントを行う回
    Returns:
        parents: 親として返された個体リスト
    """
    parents = []#親集団リスト
    for k in range(0,n):#n回のトーナメント
        i = random.randint(0,len(pop)-1)#ランダムに個体を選択
        j = random.randint(1,len(pop)-1)
        j = (i+j)%len(pop)              #もうひとつ(被らないように) 個体選択
        if pop[i].fitness > pop[j].fitness:#性能比較
                                           #※fitness(tuple)比較は最初の要素の比較になると思われる(要検証)
            parents.append(copy.deepcopy(pop[i]))#copy.deepcopy：完全なる複製の生成、
        else:                                    #               コピー元の変更の影響を受けない
            parents.append(copy.deepcopy(pop[j]))
    return parents

def binary_tournament_dom_cd(pop, n):
    """遺伝子集団のうち個体をランダムに二つ選び、性能を比較し(トーナメント)、優秀な方を親として返す
       rankで比較を行う
    Args:
        pop (gene.population): 遺伝子集団
        n (int): トーナメントを行う回
    Returns:
        parents: 親として返された個体リスト
    """
    parents = []#親集団リスト
    for k in range(0,n):#n回のトーナメント
        i = random.randint(0,len(pop)-1)#ランダムに個体を選択
        j = random.randint(1,len(pop)-1)#もうひとつ(被らないように) 個体選択
        j = (i+j)%len(pop)
        if pop[i].rank[0] != pop[j].rank[0]:#性能比較、まずrank[0]同士の比較を行い、
                                            #これが同じならrank[1]同士の比較になる
            if pop[i].rank[0] < pop[j].rank[0]:
                parents.append(copy.deepcopy(pop[i]))#copy.deepcopy：完全なる複製の生成、コピー元の変更の影響を受けない                               #               コピー元の変更の影響を受けない
            else:
                parents.append(copy.deepcopy(pop[j]))
        else: 
            if pop[i].rank[1] > pop[j].rank[1]:
                parents.append(copy.deepcopy(pop[i]))
            else:
                parents.append(copy.deepcopy(pop[j]))       
    return parents

"""交叉"""
def crossover_1p(ind1, ind2):
    """一点交叉：交差点をポイントに遺伝子を切り取りこれを入れ替える
     Args:
        ind1,ind2:親からコピーされた子供、まだただの親のクローンの個体
    Returns:
        ind1,ind2:一点交叉を行った子供個体
        nvar_change:遺伝子が何個変更されたか
    """
    k = random.randint(0,len(ind1)-1)#交叉点をランダムに選択
    xind = copy.deepcopy(ind1)
    nvar_change = 0#遺伝子が何個変更されたか
    for i in range(k,len(ind1)):
        if ind1[i] != ind2[i]:
            nvar_change += 1
        ind1[i] = ind2[i]#ind1の交差点から最後までの遺伝子をind2の変数に書き換え
        ind2[i] = xind[i]#ind2の交差点から最後までの遺伝子をind1の変数に書き換え
    return nvar_change

"""突然変異"""
def bit_flip_mutation(ind, mutp):
    """突然変異:確率で遺伝子の変数をビット反転
    Args:
        ind (gene.individual): 遺伝子個体
        mutp (float): 突然変異率(遺伝子の変数一つに対して)
    Returns:
        ind:突然変異後の遺伝子個体
        nvar_change:遺伝子が何個変更されたか
    """
    nvar_change = 0#遺伝子が何個変更されたか
    for i in range(0,len(ind)):
        if random.random() < mutp:#突然変異率に応じてビット反転
            ind[i] = (ind[i]+1)%2
            nvar_change += 1
    return nvar_change