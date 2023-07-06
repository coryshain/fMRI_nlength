import sys
import os
import gl

fROIs = [
    'all',
    'LIFGorb',
    'LIFG',
    'LMFG',
    'LAntTemp',
    'LPostTemp',
    'LAngG'
]

PDD = [
    None,
    'IFGorb',
    'IFGtri',
    None,
    'TP-aSTS',
    'pSTS',
    'TPJ'
]

PDD_all = [
    'IFGtri',
    'IFGorb',
    'TP',
    'aSTS',
    'pSTS',
    'TPJ'
]

try:
    with open('data_path.txt', 'r') as f:
        base_path = f.read().strip()
except FileNotFoundError:
    sys.stderr.write('Data path not set. Run `python -m nlength.set_data_path` before running any other scripts.\n')
    sys.stderr.flush()
    exit()

gl.resetdefaults()
gl.azimuthelevation(-90, 15)
gl.meshload('BrainMesh_ICBM152.mz3')

gl.overlayload('%s/parcels/fed_parcels.nii' % base_path)
gl.overlaycolorname(1, 'Black')
gl.overlayinvert(1, True)
gl.overlayminmax(1, 0.01, 0.01)

gl.colorbarvisible(0)
gl.overlaytransparencyonbackground(25)
gl.meshcurv()
gl.savebmp('%s/plots/evlab_parcels_all.png' % base_path)

gl.resetdefaults()
gl.azimuthelevation(-90, 15)
gl.meshload('BrainMesh_ICBM152.mz3')

gl.overlayload('%s/parcels/pdd_parcels.nii' % base_path)
gl.overlaycolorname(1, 'Black')
gl.overlayinvert(1, True)
gl.overlayminmax(1, 0.01, 0.01)

gl.colorbarvisible(0)
gl.overlaytransparencyonbackground(25)
gl.meshcurv()
gl.savebmp('%s/plots/PDD_parcels_all.png' % base_path)

for p in PDD_all:
    gl.resetdefaults()
    gl.azimuthelevation(-90, 15)
    gl.meshload('BrainMesh_ICBM152.mz3')

    gl.overlayload('%s/parcels/%s.nii' % (base_path, p))
    gl.overlaycolorname(1, 'Black')
    gl.overlayminmax(1, 0.01, 10.0)

    gl.colorbarvisible(0)
    gl.overlaytransparencyonbackground(25)
    gl.meshcurv()
    gl.savebmp('%/plots/PDD_parcels_%s.png' % (base_path, p))

for f, p in zip(fROIs, PDD):
    # LANG only
    gl.resetdefaults()
    gl.azimuthelevation(-90, 15)
    gl.meshload('BrainMesh_ICBM152.mz3')

    gl.overlayload('%s/parcels/FED_%s.nii' % (base_path, f))
    gl.overlaycolorname(1, 'Black')
    gl.overlayminmax(1, 0.01, 10.0)

    gl.colorbarvisible(0)
    gl.overlaytransparencyonbackground(25)
    gl.meshcurv()
    gl.savebmp('%s/plots/evlab_parcels_%s.png' % (base_path, f))

    # LANG-PDD overlay
    if f not in ('all', 'LMFG'):
        gl.resetdefaults()
        gl.azimuthelevation(-90, 15)
        gl.meshload('BrainMesh_ICBM152.mz3')

        gl.overlayload('%s/parcels/FED_%s.nii' % (base_path, f))
        gl.overlaycolorname(1, 'Red')

        if p is not None:
            gl.overlayload('%s/parcels/PDD_%s.nii' % (base_path, p))
            gl.overlaycolorname(2, 'Blue')

        gl.colorbarvisible(0)
        gl.overlaytransparencyonbackground(25)
        gl.meshcurv()
        gl.savebmp('%s/plots/evlab_pdd_parcels_%s.png' % (base_path, f))

