from pathlib import Path
from tqdm import tqdm
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas
import numpy as np

import pybdsim

def _generate_test_outputs(dirToScan) :
    root = Path(dirToScan)
    files = list(root.rglob("*.root"))
    return files

def generate(dirToScan, nFileLimit = 10000):
    files = _generate_test_outputs(dirToScan)
    filesSkipped = []
    filesNoRun = []

    dd = {}
    dd['index'] = []
    dd['file_path'] = []
    dd['file_name'] = []
    dd['nGenerate'] = []
    dd['memoryUsageMb'] = []
    dd['memoryUsageMbStd'] = []
    dd['memoryUsageMbMin'] = []
    dd['memoryUsageMbMax'] = []

    dd['eventDurationCPU'] = []
    dd['eventDurationCPUStd'] = []
    dd['eventDurationCPUMin'] = []
    dd['eventDurationCPUMax'] = []

    dd['eventDurationWall'] = []
    dd['eventDurationWallStd'] = []
    dd['eventDurationWallMin'] = []
    dd['eventDurationWallMax'] = []

    dd['runDurationCPU'] = []
    dd['runDurationWall'] = []

    for i,f in tqdm(enumerate(files)) :

        if i == nFileLimit :
            break

        if "fodo-no-primaries.root" in str(f) or \
           "shower_cyl_sph.root" in str(f) :
            continue

        #print(i,f)
        #print(pybdsim.DataPandas._ROOTFileType(str(f)))

        p = pybdsim.DataPandas.Load(str(f))

        if not p :
            continue

        # options
        o = p.get_options()

        # run
        r = p.get_run()

        # event summary
        es = p.get_event_summary()

        # output data
        dd['file_path'].append(str(f))
        dd['file_name'].append(str(Path(f.name).stem))

        # index
        dd['index'].append(i)

        # number of events
        dd['nGenerate'].append(o['nGenerate'][0])

        # memory usage
        dd['memoryUsageMb'].append(es['memoryUsageMb'].mean())
        dd['memoryUsageMbStd'].append(es['memoryUsageMb'].std())
        dd['memoryUsageMbMin'].append(es['memoryUsageMb'].min())
        dd['memoryUsageMbMax'].append(es['memoryUsageMb'].max())

        # duration CPU
        dd['eventDurationCPU'].append(es['durationCPU'].mean())
        dd['eventDurationCPUStd'].append(es['durationCPU'].std())
        dd['eventDurationCPUMin'].append(es['durationCPU'].min())
        dd['eventDurationCPUMax'].append(es['durationCPU'].max())

        # duration Wall
        dd['eventDurationWall'].append(es['durationWall'].mean())
        dd['eventDurationWallStd'].append(es['durationWall'].std())
        dd['eventDurationWallMin'].append(es['durationWall'].min())
        dd['eventDurationWallMax'].append(es['durationWall'].max())

        # run duration CPU
        try :
            dd['runDurationCPU'].append(r['durationCPU'][0])
        except :
            print(str(f))
            dd['runDurationCPU'].append(0)

        # run duration Wall
        try :
            dd['runDurationWall'].append(r['durationWall'][0])
        except :
            print(str(f))
            dd['runDurationWall'].append(0)

    return pandas.DataFrame(dd)

def plot_all(df) :
    ranges = list(range(0,len(df), 50))
    ranges.append(len(df))
    print(ranges)

    for i,r in enumerate(ranges[0:-1]) :
        print(i, ranges[i], ranges[i+1])
        plot_page(df, i, ranges[i], ranges[i+1])


    with PdfPages("plot_all.pdf") as pdf:
        for fig_num in plt.get_fignums():
            pdf.savefig(plt.figure(fig_num))


def plot_selection(df, pathContains ="", fileContains = "") :
    pass

def plot_page(df, ifigure=1, istart = 0, iend=50) :
    a4size = np.array([8.27, 11.69])*0.7

    plt.figure(figsize=tuple(a4size)) # A4

    plt.subplot(1,4,1)
    plt.barh(df['file_name'][istart:iend], df['nGenerate'][istart:iend])
    # plt.xlim(0, df['nGenerate'].max())

    plt.subplot(1,4,2)
    plt.barh(df['file_name'][istart:iend].str.pad(width=30, side='right'),
             df['runDurationWall'][istart:iend])
    plt.gca().set_yticks([])
    plt.xlim(0, df['runDurationWall'].max())

    plt.subplot(1,4,3)
    plt.barh(df['file_name'][istart:iend], df['memoryUsageMb'][istart:iend])
    plt.gca().set_yticks([])
    plt.xlim(0, df['memoryUsageMb'].max())


    plt.tight_layout()