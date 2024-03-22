import re
import os
from collections import defaultdict


def tokenize(fname, pattern_type=0):
    if pattern_type == 0:
        pattern = 'stats_Run([.\d]+)_(.+)_([.\d]+)(:?V|v)_([.\d]+)(:?keV|KeV|kev)_(.+)_([.\d]+).root'
    elif pattern_type == 1:
        pattern = 'stats_Run([.\d]+)_(.+)_([.\d]+)(:?V|v)_([.\d]+)(:?keV|KeV|kev)_([.\d]+).root'
    elif pattern_type == 2:
        pattern = 'stats_Run([.\d]+)_(.+)_([.\d]+)(:?V|v)_(.+)([.\d]+)(:?keV|KeV|kev)_([.\d]+).root'
    elif pattern_type == 3:
        pattern = 'stats_Run([.\d]+)_(.+)_([.\d]+)(:?V|v)_([0-9]+)(:?keV|KeV|kev).root'
    elif pattern_type == 4:
        pattern = (
            'stats_Run([.\d]+)_(.+)_([.\d]+)(:?V|v)_([0-9]+)(:?keV|KeV|kev)_(.+).root'
        )
    else:
        return None
    c_pattern = re.compile(pattern)

    tokens = c_pattern.match(fname)

    if tokens is None:
        return tokenize(fname, pattern_type + 1)

    return tokens


def resolve_filename(files, pattern_type=0):
    # pattern_0 = r"stats_Run([.\d]+)_(.+)_([.\d]+)(:?V|v)_([.\d]+)(:?keV|KeV|kev)_(.+?:?)_(:?ch)"
    # pattern_1 = r"stats_Run([.\d]+)_(.+)_([.\d]+)(:?V|v)_([.\d]+)(:?keV|KeV|kev)_(.+?)(:?[.\d]+)_(:?ch)"

    group_map = {
        "run": 1,
        "sensor": 2,
        "voltage": 3,
        "volt_unit": 4,
        "energy": 5,
        "energy_unit": 6,
        "user_note": 7,
    }

    filesdict = defaultdict(list)
    infodict = {}
    for f in files:
        tokens = tokenize(os.path.basename(f))
        if tokens is None:
            # raise ValueError(f"unable to parse {f}")
            print(f"unable to parse {f}")
            continue
        run = tokens.group(group_map['run'])
        filesdict[run].append(f)
        if run not in infodict:
            temp = {}
            for k, v in group_map.items():
                try:
                    temp[k] = tokens.group(v)
                except:
                    print(f"cannot find {k}")
                    temp[k] = None
            infodict[run] = temp

    return filesdict, infodict
