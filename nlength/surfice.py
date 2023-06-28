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

gl.resetdefaults()
gl.azimuthelevation(-90, 15)
gl.meshload('BrainMesh_ICBM152.mz3')

gl.overlayload('../parcels/fed_parcels.nii')
gl.overlaycolorname(1, 'Black')
gl.overlayinvert(1, True)
gl.overlayminmax(1, 0.01, 0.01)

gl.colorbarvisible(0)
gl.overlaytransparencyonbackground(25)
gl.meshcurv()
gl.savebmp('../plots/evlab_parcels_all.png')

gl.resetdefaults()
gl.azimuthelevation(-90, 15)
gl.meshload('BrainMesh_ICBM152.mz3')

gl.overlayload('../parcels/pdd_parcels.nii')
gl.overlaycolorname(1, 'Black')
gl.overlayinvert(1, True)
gl.overlayminmax(1, 0.01, 0.01)

gl.colorbarvisible(0)
gl.overlaytransparencyonbackground(25)
gl.meshcurv()
gl.savebmp('../plots/PDD_parcels_all.png')

for p in PDD_all:
    gl.resetdefaults()
    gl.azimuthelevation(-90, 15)
    gl.meshload('BrainMesh_ICBM152.mz3')

    gl.overlayload('../parcels/%s.nii' % p)
    gl.overlaycolorname(1, 'Black')
    gl.overlayminmax(1, 0.01, 10.0)

    gl.colorbarvisible(0)
    gl.overlaytransparencyonbackground(25)
    gl.meshcurv()
    gl.savebmp('../plots/PDD_parcels_%s.png' % p)

for f, p in zip(fROIs, PDD):
    # LANG only
    gl.resetdefaults()
    gl.azimuthelevation(-90, 15)
    gl.meshload('BrainMesh_ICBM152.mz3')

    gl.overlayload('../parcels/FED_%s.nii' % f)
    gl.overlaycolorname(1, 'Black')
    gl.overlayminmax(1, 0.01, 10.0)

    gl.colorbarvisible(0)
    gl.overlaytransparencyonbackground(25)
    gl.meshcurv()
    gl.savebmp('../plots/evlab_parcels_%s.png' % f)

    # LANG-PDD overlay
    if f not in ('all', 'LMFG'):
        gl.resetdefaults()
        gl.azimuthelevation(-90, 15)
        gl.meshload('BrainMesh_ICBM152.mz3')

        gl.overlayload('../parcels/FED_%s.nii' % f)
        gl.overlaycolorname(1, 'Red')

        if p is not None:
            gl.overlayload('../parcels/PDD_%s.nii' % p)
            gl.overlaycolorname(2, 'Blue')

        gl.colorbarvisible(0)
        gl.overlaytransparencyonbackground(25)
        gl.meshcurv()
        gl.savebmp('../plots/evlab_pdd_parcels_%s.png' % f)

