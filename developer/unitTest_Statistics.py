from pathlib import Path
from tqdm import tqdm
import pandas

import pybdsim

def _generate_test_outputs(dirToScan) :
    root = Path(dirToScan)
    files = list(root.rglob("*.root"))
    return files

def generate(dirToScan):
    files = _generate_test_outputs(dirToScan)
    filesSkipped = []

    dd = {}
    dd['index'] = []
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
        dd['file_name'].append(str(f))

        # index
        dd['index'].append(i)

        # number of events
        dd['nGenerate'].append(o['nGenerate'])

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
        dd['runDurationCPU'].append(r['durationCPU'][0])

        # run duration Wall
        dd['runDurationWall'].append(r['durationWall'][0])


    return pandas.DataFrame(dd)

def plot(df) :
    pass

