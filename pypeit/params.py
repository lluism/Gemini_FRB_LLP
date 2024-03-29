from astropy.coordinates import SkyCoord
from astropy import units as u

inst_map = {'gmos_south': 'gemini_gmos_south_ham',
            'gmos_north': 'gemini_gmos_north_ham'}

params = {'gmos_south': ['[calibrations]\n',
'    [[biasframe]]\n',
'        [[[process]]]\n',
'            combine = median\n',
'            use_biasimage = False\n',
'            shot_noise = False\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[darkframe]]\n',
'        [[[process]]]\n',
'            mask_cr = True\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[arcframe]]\n',
'        [[[process]]]\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[tiltframe]]\n',
'        [[[process]]]\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[pixelflatframe]]\n',
'        [[[process]]]\n',
'            combine = median\n',
'            satpix = nothing\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[alignframe]]\n',
'        [[[process]]]\n',
'            satpix = nothing\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[traceframe]]\n',
'        [[[process]]]\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[illumflatframe]]\n',
'        [[[process]]]\n',
'            satpix = nothing\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[skyframe]]\n',
'        [[[process]]]\n',
'            mask_cr = True\n',
'            lamaxiter = 3\n',
'            sigclip = 3.0\n',
'            noise_floor = 0.01\n',
'    [[standardframe]]\n',
'        [[[process]]]\n',
'            mask_cr = True\n',
'            lamaxiter = 3\n',
'            sigclip = 3.0\n',
'            noise_floor = 0.01\n',
'    [[wavelengths]]\n',
'        method = full_template\n',
'        lamps = CuI, ArI, ArII\n',
'        rms_threshold = 0.4\n',
'        nsnippet = 1\n',
'    [[slitedges]]\n',
'        fit_order = 3\n',
'        bound_detector = True\n',
'    [[tilts]]\n',
'        tracethresh = 10.0\n',
'[scienceframe]\n',
'    [[process]]\n',
'        mask_cr = True\n',
'        noise_floor = 0.01\n',
'[flexure]\n',
'    spec_method = boxcar\n',
'[sensfunc]\n',
'    multi_spec_det = 1, 2, 3\n',
'    algorithm = IR\n',
'    [[IR]]\n',
'        telgridfile = TelFit_LasCampanas_3100_26100_R20000.fits\n'],
'gmos_north':  ['[calibrations]\n',
'    [[biasframe]]\n',
'        [[[process]]]\n',
'            combine = median\n',
'            use_biasimage = False\n',
'            shot_noise = False\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[darkframe]]\n',
'        [[[process]]]\n',
'            mask_cr = True\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[arcframe]]\n',
'        [[[process]]]\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[tiltframe]]\n',
'        [[[process]]]\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[pixelflatframe]]\n',
'        [[[process]]]\n',
'            combine = median\n',
'            satpix = nothing\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[alignframe]]\n',
'        [[[process]]]\n',
'            satpix = nothing\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[traceframe]]\n',
'        [[[process]]]\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[illumflatframe]]\n',
'        [[[process]]]\n',
'            satpix = nothing\n',
'            use_pixelflat = False\n',
'            use_illumflat = False\n',
'    [[skyframe]]\n',
'        [[[process]]]\n',
'            mask_cr = True\n',
'            lamaxiter = 3\n',
'            sigclip = 3.0\n',
'            noise_floor = 0.01\n',
'    [[standardframe]]\n',
'        [[[process]]]\n',
'            mask_cr = True\n',
'            lamaxiter = 3\n',
'            sigclip = 3.0\n',
'            noise_floor = 0.01\n',
'    [[wavelengths]]\n',
'        method = full_template\n',
'        lamps = CuI, ArI, ArII\n',
'        rms_threshold = 0.4\n',
'        nsnippet = 1\n',
'    [[slitedges]]\n',
'        fit_order = 3\n',
'        bound_detector = True\n',
'    [[tilts]]\n',
'        tracethresh = 10.0\n',
'[scienceframe]\n',
'    [[process]]\n',
'        mask_cr = True\n',
'        noise_floor = 0.01\n',
'[flexure]\n',
'    spec_method = boxcar\n',
'[sensfunc]\n',
'    multi_spec_det = 1, 2, 3\n',
'    algorithm = IR\n',
'    [[IR]]\n',
'        telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits\n'],}



class StandardStar():
    def __init__(self,ra='',dec=''):
        self.coord = SkyCoord('{} {}'.format(ra,dec),frame='icrs',
            unit=(u.hourangle, u.deg))

def construct_standard_star_library():
    ''' Library of standard stars '''

    ssl = {
        'bd174708': StandardStar(ra='22:11:31.38', dec='+18:05:34.2'),
        'bd262606': StandardStar(ra='14:49:02.36', dec='+25:42:09.1'),
        'bd284211': StandardStar(ra='21:51:11.02', dec='+28:51:50.4'),
        'bd332642': StandardStar(ra='15:51:59.89', dec='+32:56:54.3'),
        'feige15': StandardStar(ra='01:49:09.49', dec='+13:33:11.8'),
        'feige24': StandardStar(ra='02:35:07.59', dec='+03:43:56.8'),
        'feige25': StandardStar(ra='02:38:37.79', dec='+05:28:11.3'),
        'feige34': StandardStar(ra='10:39:36.74', dec='+43:06:09.2'),
        'feige56': StandardStar(ra='12:06:47.24', dec='+11:40:12.7'),
        'feige92': StandardStar(ra='14:11:31.88', dec='+50:07:04.1'),
        'feige98': StandardStar(ra='14:38:15.75', dec='+27:29:32.9'),
        'feige110':StandardStar(ra='23:19:58.4', dec='-05:09:56.2'),
        'g158100': StandardStar(ra='00:33:54', dec='-12:07:57'), ###
        'g191b2b': StandardStar(ra='05:05:30.62', dec='+52:49:51.9'),
        'gd71': StandardStar(ra='05:52:27.62', dec='+15:53:13.2'),
        'gd248': StandardStar(ra='23:26:07', dec='+16:00:21'), ###
        'hd19445': StandardStar(ra='03:08:25.59', dec='+26:19:51.4'),
        'hd84937':StandardStar(ra='09:48:56.1',dec='+13:44:39.3'),
        'hz43':  StandardStar(ra='13:16:21.85', dec='+29:05:55.4'),
        'hz44': StandardStar(ra='13:23:35.26', dec='+36:07:59.5'),
        'ltt0379':StandardStar(ra='18:36:25.95', dec='-44:18:36.9'),
        'ltt1020':StandardStar(ra='01:54:50.27',dec='-27:28:35.7'),
        'ltt1788': StandardStar(ra='03:48:22.61', dec='-39:08:37.0'),
        'ltt2415': StandardStar(ra='05:56:24.74', dec='-27:51:32.4'),
        'ltt3218': StandardStar(ra='08:41:32.43', dec='-32:56:32.9'),
        'ltt3864': StandardStar(ra='10:32:13.62', dec='-35:37:41.7'),
        'ltt4364': StandardStar(ra='11:45:42.92', dec='-64:50:29.5'),
        'ltt6248': StandardStar(ra='15:38:59.648', dec='-28:35:36.97'),
        'ltt7379': StandardStar(ra='18:36:25.950', dec='-44:18:36.90'),
        }
    return ssl
