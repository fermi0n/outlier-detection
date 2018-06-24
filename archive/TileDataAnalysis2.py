import numpy as np
import scipy.stats as ss
import scipy.spatial.distance as ssd
from core import kshape, zscore, _sbd
from operator import itemgetter
import DataAnalysis as da
import CalsLoader
from TileDataAnalysis import DataAnalyserByTile
#Uncomment below line to use DTW
#import mlpy as mlpy
import matplotlib.pyplot as plt
import argparse
import pickle


if __name__ == '__main__':

    #Get plotting ready
    colours = ['#AE70ED','#FFB60B','#62A9FF','#59DF00']

    sp = 0
    ppl = 16

    maxv = 10.0
    #    maxv = max(max([max(tiledata.allx_obs_dict[obs][ii]) for ii in range(128)]), max([max(tiledata.ally_obs_dict[obs][ii]) for ii in range(128)]))

    fig = plt.figure(figsize=(18.0, 10.0))


    #Unpickle this data
    tile = '4'

    with open ('/lustre/projects/p048_astro/dmendonca/%s_data.pickle'%(tile), 'rb') as f:
        tileloader = pickle.load(f)

    #distancedata = (hardcoded distances)
    distancedata = [['1127245592', 4, 9.5602952740780776], ['1135431024', 4, 5.98450708966379], ['1135431144', 4, 5.9819937072337144], ['1132846192', 4, 5.7539484278182513], ['1135431272', 4, 5.6430759584731582], ['1128365600', 4, 5.5488758604348885], ['1135430904', 4, 5.5275604172468942], ['1129399688', 4, 5.4831804882871449], ['1128366576', 4, 5.4475486704694438], ['1132847168', 4, 5.4397316320332516], ['1129400544', 4, 5.4187812917066021], ['1132846072', 4, 5.4187062574464475], ['1128365480', 4, 5.414292628174751], ['1132846312', 4, 5.4033280771055585], ['1129399568', 4, 5.3854558923712146], ['1132847048', 4, 5.3693085069402242], ['1128365720', 4, 5.3624400039987687], ['1129401032', 4, 5.3513340060876455], ['1125953008', 4, 5.3464655937679844], ['1128194248', 4, 5.3243430749864178], ['1132847536', 4, 5.316100614581651], ['1125953128', 4, 5.2896534168313529], ['1132847416', 4, 5.2586894778321529], ['1128366696', 4, 5.2581366689605957], ['1129399448', 4, 5.2509852480326611], ['1125953984', 4, 5.2397407303320822], ['1125953248', 4, 5.2370265517456431], ['1132847288', 4, 5.2184808193957641], ['1126126312', 4, 5.2116512566051139], ['1128366944', 4, 5.2110826102203553], ['1128194736', 4, 5.2045348251180101], ['1126298640', 4, 5.1864865075141529], ['1128366816', 4, 5.1852085137230768], ['1129399816', 4, 5.1812145496059507], ['1125952888', 4, 5.1775475673044422], ['1125954472', 4, 5.1561546720663447], ['1132845952', 4, 5.1535954563527984], ['1135431880', 4, 5.151103224129181], ['1126125336', 4, 5.1454622122568079], ['1126126800', 4, 5.0961661235539193], ['1128365360', 4, 5.0896400890167799], ['1126125456', 4, 5.0746542788807147], ['1126125576', 4, 5.068177003298616], ['1129400424', 4, 5.0508415031261897], ['1132847656', 4, 5.0443386754195121], ['1128365848', 4, 5.036268348903695], ['1129400912', 4, 5.0345384727055764], ['1128194496', 4, 5.0159731247934154], ['1129400792', 4, 5.0151131991590985], ['1126299128', 4, 4.9977275654346709], ['1129399328', 4, 4.9926265026508201], ['1126125216', 4, 4.9839128576067377], ['1129400664', 4, 4.9808918811963201], ['1125954104', 4, 4.9697743134197463], ['1128366456', 4, 4.9613888800177559], ['1128367184', 4, 4.9519383913233899], ['1126126192', 4, 4.9508269648677814], ['1132846928', 4, 4.9492480441173043], ['1125954344', 4, 4.9182294897935162], ['1135431392', 4, 4.9082391257658289], ['1135431632', 4, 4.9075951558625537], ['1128194368', 4, 4.8975677024095896], ['1128194616', 4, 4.8866800205850112], ['1125954224', 4, 4.8302708748828964], ['1125953856', 4, 4.8301759815778658], ['1135431512', 4, 4.8133404916718376], ['1135431760', 4, 4.8076905148988267], ['1126298880', 4, 4.8012426916069417], ['1128194128', 4, 4.8003469752714505], ['1126126680', 4, 4.8001215748858153], ['1126298760', 4, 4.7959041995217895], ['1126126552', 4, 4.7779845101178244], ['1126299008', 4, 4.7755561772589221], ['1129401152', 4, 4.769909358259584], ['1126298520', 4, 4.758368363764454], ['1132846440', 4, 4.7479272367560563], ['1126126432', 4, 4.7177656623194686], ['1132846680', 4, 4.6794717806822206], ['1132846800', 4, 4.6599764723819099], ['1128366208', 4, 4.6592412793755082], ['1129400304', 4, 4.6331940458737542], ['1129399936', 4, 4.6255131860149943], ['1128366336', 4, 4.6120095555476972], ['1128194856', 4, 4.5752051019363789], ['1129400176', 4, 4.5656473760505394], ['1127245472', 4, 4.5654116492032371], ['1125954592', 4, 4.5519999192118199], ['1126299248', 4, 4.5460555850917546], ['1127245712', 4, 4.5425610889345105], ['1127245344', 4, 4.5416559569309811], ['1128193880', 4, 4.5055048507081068], ['1125953616', 4, 4.4959008421965159], ['1125953736', 4, 4.464804863084245], ['1128366088', 4, 4.4644687684326128], ['1126125944', 4, 4.4597866319139943], ['1128365968', 4, 4.456449014887693], ['1129400056', 4, 4.443335756591595], ['1126126064', 4, 4.4296663241279477], ['1125953376', 4, 4.4214273066483836], ['1126126920', 4, 4.4071410734228627], ['1128193760', 4, 4.4008123857394983], ['1128194000', 4, 4.3906314084931104], ['1126125704', 4, 4.3493729492391076], ['1126125824', 4, 4.3460066893953231], ['1125953496', 4, 4.3459274567659705], ['1126298392', 4, 4.3121738257101896], ['1127245832', 4, 4.311435345527002], ['1132846560', 4, 4.3088127208339637], ['1126298152', 4, 4.2636351864176962], ['1126211376', 4, 4.2409417951286983], ['1126298272', 4, 4.2298547540752205], ['1126211256', 4, 3.9259378098575413], ['1127245224', 4, 3.679425054747437], ['1127245952', 4, 3.6260979581757655], ['1129487320', 4, 3.5842587031467703], ['1129485976', 4, 3.3851901703610148], ['1129487200', 4, 3.2680975870415225], ['1126039416', 4, 3.2669865276672212], ['1126211744', 4, 3.2629306018215964], ['1129486584', 4, 3.2495618968989604], ['1129486712', 4, 3.227837523124057], ['1129486832', 4, 3.1529775043749964], ['1129485856', 4, 3.0834150819079378], ['1129486464', 4, 3.0580280010791085], ['1129486952', 4, 3.0311686627500443], ['1135173016', 4, 2.997864307204352], ['1126213088', 4, 2.9769607434375329], ['1124316008', 4, 2.9542537405426423], ['1129487072', 4, 2.9018473228697421], ['1129486096', 4, 2.8955823298198857], ['1126212472', 4, 2.8887069285057025], ['1129486344', 4, 2.8837352030165029], ['1126211984', 4, 2.8385879616259642], ['1129486224', 4, 2.8017141667439418], ['1132070840', 4, 2.7495563862355743], ['1132759912', 4, 2.7342635233703931], ['1135346448', 4, 2.7298472845374353], ['1126211624', 4, 2.7248146010240655], ['1127247056', 4, 2.6810996681719175], ['1129314136', 4, 2.6613561711546838], ['1126212112', 4, 2.6483232549424627], ['1105360912', 4, 2.6452015218262996], ['1127246320', 4, 2.6339497898909552], ['1126212352', 4, 2.6282694837381895], ['1126039656', 4, 2.6178786630633355], ['1105016376', 4, 2.5893768057823734], ['1126040024', 4, 2.5888031393808535], ['1125006792', 4, 2.5836610916746254], ['1124316744', 4, 2.5817666435746931], ['1124316864', 4, 2.5786621360991329], ['1124144168', 4, 2.5682876150174776], ['1126040752', 4, 2.5636789108875586], ['1124316136', 4, 2.5556670754843056], ['1129314872', 4, 2.554051524625685], ['1105016016', 4, 2.5430428170216839], ['1129485736', 4, 2.537961479601778], ['1124143808', 4, 2.5364686067300504], ['1126040144', 4, 2.5364613455694878], ['1124920624', 4, 2.5348957325360124], ['1105189072', 4, 2.5305608840853444], ['1105447568', 4, 2.5298693704151609], ['1125006056', 4, 2.5292885336185349], ['1129314256', 4, 2.5241579663175373], ['1126212600', 4, 2.5211540915228077], ['1124919888', 4, 2.5137276051692616], ['1124143680', 4, 2.5110283965451994], ['1126212960', 4, 2.5104206417350081], ['1126212720', 4, 2.5051735735516156], ['1129314016', 4, 2.5010563726635695], ['1124316256', 4, 2.4993825453424812], ['1105361400', 4, 2.4978534002584873], ['1129314992', 4, 2.4970359966407703], ['1125005320', 4, 2.4950116799399931], ['1127246440', 4, 2.4945204411378783], ['1126039296', 4, 2.4921247530430755], ['1105619160', 4, 2.4911501088587462], ['1125006304', 4, 2.490301855320955], ['1126211864', 4, 2.4900738993182352], ['1129313528', 4, 2.4867547970534871], ['1126040264', 4, 2.4852529107319832], ['1126211496', 4, 2.4846257264617906], ['1124230576', 4, 2.4843944764534966], ['1127246200', 4, 2.4828072008743787], ['1105361160', 3, 3.296515069293386], ['1127246688', 3, 3.2837720013470655], ['1105360912', 3, 3.2683498786802669], ['1126040392', 3, 3.2577381822823916], ['1126212600', 3, 3.2497315853780346], ['1124316864', 3, 3.2472876796524881], ['1125006424', 3, 3.2470478281034638], ['1124920256', 3, 3.2439344577461915]]

    for obsid, tile, distance in distancedata[:127]:

        ax = fig.add_subplot(8,16,sp+1,)
        ax.plot(tileloader.allx_obs_dict[obsid], colours[1])
        ax.plot(tileloader.ally_obs_dict[obsid], colours[2])
        plt.title('Obs %s\n, distance %.1f'%(obsid, distance), fontsize=6)

        if ppl != 16:
            plt.setp(ax.get_xticklabels(), visible=False) # plot setup
            plt.setp(ax.get_yticklabels(), visible=False)

        if ppl == 16:
            ppl = 0
            plt.setp(ax.get_xticklabels(), visible=False)

        ppl += 1

        plt.ylim([-0.1,maxv])

        if sp == 15:
            XX_amp, = ax.plot([],[], colours[1],label='XX',linewidth=3.0)
            YY_amp, = ax.plot([],[], colours[2],label='YY',linewidth=3.0)
            ax.legend((XX_amp,YY_amp),('XX','YY'), bbox_to_anchor=(0, 2, .12, .12),prop={'size':14})

        sp += 1

        plt.tight_layout()
        fig.subplots_adjust(top=0.9)
        plt.suptitle('Amps | Tile %s' %(tile),fontsize=18)

    plt.show()