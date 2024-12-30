import random, copy
from . import ea_base as ea
from . import moea_base as moea
from .tool.fitness_method_match import match_ind_pop as mip
import os
from .tool.create_dtinfo import store_DTinfo

def nsgaii_survival_selection(pop, popsize):
    """fitnessから非支配レベルfrontと混雑距離cdを計算しソーティング、生存選択を行う
    Args:
        pop (population): 遺伝子集団
        popsize (int): 生き残らせる個体数
    """
    
    nobj = pop[0].nobj#評価関数の数
    
    pop_inf=ea.Population()###実行不可能解を別処理とする
    size=len(pop)
    for i in range(size-1, -1, -1):
        if pop[i].penalty>0:
            pop_inf.insert(0, pop[i])#解を挿入
            pop.remove(pop[i])#元の集団から解は消え去る
    
    if len(pop)!=0:
        fronts = moea.non_dominated_sorting(pop)#fronts:list　0番目にはflont0の個体集団、1番目にはflont1の個体集団...が入ってる
    else:
        fronts=ea.Population()
    #このときpopは空
    
    
    #popsizeのキャパを超えるまでfront0の個体集団,front1の個体集団...という順に空のpopに同一個体集団まるまる詰め込む
    i = 0
    if len(fronts)!=0:
        while len(pop) < popsize and len(pop) + len(fronts[0]) <= popsize:
            fi = fronts.pop(0)#frontsの0番目の要素を返し、さらにこれをリストから消去
            moea.front_rank(fi, i)#fi(同フロント集団)の個体のindividual.rank[0]にフロントランクを書き込む
            moea.crowding_distance(fi, nobj)#fiの個体のindividual.rank[1]に混雑距離を書き込む
            pop.extend(fi)#pop(list)の末尾にfi(list)をくっつける
            i += 1#フロントが一つ下がることを意味する
            if len(fronts)==0:
                break
        
    
    #popsizeのキャパを超えたら同一フロント内で混雑距離でソーティングし足きりを行う
    if len(pop) < popsize and len(fronts)!=0:
        fi = fronts.pop(0)
        moea.front_rank(fi, i)
        moea.crowding_distance(fi, nobj)
        fi.sort(key = lambda ind: ind.rank[1], reverse = True)
        pop.extend(fi[: popsize - len(pop)])
    
    #実行不可能解の処理、フロントを同一で１００００とする、またCDを書き込むところをペナルティを負に変換した値を書き込む
    if len(pop) < popsize:
        moea.front_rank(pop_inf,10000)
        for ind in pop_inf:
            ind.rank[1]=-ind.penalty
        pop_inf.sort(key = lambda ind: ind.rank[1],reverse = True)
        pop.extend(pop_inf[: popsize - len(pop)])
        

    if len(pop) != popsize:#エラー検知
        print("Somethig is wrong, len(pop) = " + 
              str(len(pop)) + ", popsize = " + str(popsize))
    
def nsgaii(evaluate=None, select=None, recombine=None, mutate=None, 
    seed=None, psize=None, nobj=None, nvar=None, vlow=None, vhigh=None, 
    initype=None, ngen=None, pcx=None, pmut=None, keepclones = False):
    """NSGA-llのアルゴリズムに基づいた処理を行う
    Args:
        evaluate(function):遺伝子評価の関数 
        select(function):親選択の関数
        recombine(function):交叉の関数
        mutate(function):突然変異の関数
        seed(int):遺伝子集団のランダムな初期集団生成時(及び交叉確率)のシード値
        psize(int):遺伝子集団の個体の数
        nvar(int):遺伝子の変数の数
        vlow,vhigh(float):遺伝子の変数の範囲(vlow~vhigh)、initypeがbinaryならvhighは初期個体生成において一つの変数が1になる確率
        initype(str):変数のタイプ
        ngen(int):世代数
        pcx(float):交叉確率
        pmut(float):突然変異確率(一つの変数に対して)
        keepclones(bool):親と全く同じ遺伝子を残すか否か
    """
    random.seed(seed)
    """ Initial population  """#初期集団生成
    pop = ea.Population(size=psize, nobj=nobj, nvar=nvar, vlow=vlow, 
                        vhigh=vhigh, initype=initype)
    """ Evaluate the initial population """
    for ind in pop:
        ###変更
        fitness,columns,values,dtinfo,mt,penalty= evaluate(ind)
        ###目的関数、csvに書き込むカラムリスト、それに対応する値リスト、決定木情報、遺伝子一致率を評価関数に採用するか否か
        if mt==1:
            match=mip(ind,pop)
            fitness.append(-match)#遺伝子一致率
            values.append(match)
        else:
            values.append(0)
        ind.fitness=tuple(fitness)
        ind.values=values#書き込む値を個体に保持
        ind.dtinfo=dtinfo#遺伝子からつくられる決定木の内部構造情報を保持
        ind.penalty=penalty#ペナルティを保存０より大きければ実行不可能解

        
    """ Output the population """#初期集団output
    nsgaii_survival_selection(pop, psize)#評価によるソーティング(非支配レベルと混雑距離）+下位個体を淘汰？
    pop.fprintpopcsv("pop_g0.csv",columns)#csvファイルとしてpddataframeで書き込む
    os.mkdir("clf_info_gen0")#決定木情報を保存する
    for i in range(len(pop)):#
        store_DTinfo(pop[i].dtinfo,"clf_info_gen0/clf"+str(i)+".txt")

    
    print(' --Generation 0')
   
    for g in range(1, ngen+1):
        """ Select the next generation individuals """#親選択　
        parents = select(pop, psize)#ea_base参照(親選択法)
        """ Clone the selected individuals """#親のクローンを生成
        offspring = copy.deepcopy(parents)

        """ Apply crossover and mutation on the offspring """#クローン個体を二つ選んで交叉 突然変異
        for ind1, ind2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < pcx:
                recombine(ind1, ind2)#ea_base参照(交叉法)
            mutate(ind1, pmut)#ea_base参照(突然変異法)
            mutate(ind2, pmut)

        """ Evaluate offspring population """#子集団(クローン個体)を評価
        """変更"""
        for ind in offspring:
            fitness,columns,values,dtinfo,mt,penalty = evaluate(ind)
            if mt==1:
                match=mip(ind,pop+offspring)
                fitness.append(-match)#遺伝子一致率
                values.append(round(match,3))
            else:
                values.append(0)
            ind.fitness=tuple(fitness)
            ind.values=values
            ind.dtinfo=dtinfo
            ind.penalty=penalty
            
        """ Delete clones """#親と全く同じ遺伝子は抹殺
        if keepclones == False:
            join_pop = []
            for ind in (pop+offspring):
                if ind not in join_pop:
                    join_pop.append(ind)
            pop[:] = join_pop[:]
        else:
            pop.extend(offspring)
        
        """ Truncate the population (extinctive selection) 
            The new population is the best among pop and offspring
        """

        nsgaii_survival_selection(pop, psize)#評価ソート＆淘汰
      
        """ Output the population """
        if g%(ngen/10) == 0:
            title="pop_g" + str(g) + ".csv"
            pop.fprintpopcsv(title,columns)
        
        if g%(ngen/5) == 0:
            os.mkdir("clf_info_gen"+str(g))
            for i in range(len(pop)):
                store_DTinfo(pop[i].dtinfo,"clf_info_gen"+str(g)+"/clf"+str(i)+".txt")
                        
        if g%(ngen/4) == 0:
            print(' --Generation ', g)
        
    print('Ends Nsgaii')
    return pop