#
# Copyright (c) European Synchrotron Radiation Facility (ESRF)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

from __future__ import division, print_function
import os
import sys
import json
import jsonschema
import pathlib
import numpy as np
import scipy.spatial as sp

from cctbx import miller
from cctbx import crystal

from edna2.tasks.AbstractTask import AbstractTask
from edna2.utils import UtilsLogging

__author__ = ['S. Basu']
__license__ = 'MIT'
__date__ = '20/11/2020'

"Inspired from Computing workshop in ECM30 conference: https://github.com/keitaroyam/ECACOMSIG-ECM30"

logger = UtilsLogging.getLogger()


class ExecSpotPrediction(AbstractTask):

    def __init__(self, jData):
        AbstractTask.__init__(self, inData=jData)
        self.xparamdict = self.digestxparamFile(inData=jData)

    def getInDataSchema(self):
        return {
            "image": {"type": "string"},
            "XDSspotFile": {"type", "string"},
            "XDSidxrefFile": {"type", "string"},
            "XDSxparamFile": {"type", "string"},
            "resolution": {"type", "string"},
            "frame": {"type", "integer"},
            "mosaicity": {"type", "number"}
        }

    @staticmethod
    def nearest_neighbor(spot_array, predicted_array, dist_limit, ofh):
        kdtree = sp.KDTree(spot_array)
        nNeighbor = 0
        for (xcal, ycal, phi, zeta) in predicted_array:
            d, idx = kdtree.query((xcal, ycal), k=1, p=2)
            print(d, idx)
            if idx == kdtree.n:
                continue
            if d < dist_limit:
                ofh.write("%.2f %.2f \n" % tuple(spot_array[idx]))
                nNeighbor += 1

        return nNeighbor

    @staticmethod
    def digestxparamFile(inData):
        xparaminfo = dict()
        idxfile = pathlib.Path(inData['XDSxparamFile'])
        if not idxfile.exists():
            return xparaminfo
        idx = open(idxfile, 'r')
        allLines = idx.readlines()
        for ii in range(len(allLines)):
            if ii == 1:
                line = allLines[ii].split()
                xparaminfo['frame_no'] = int(line[0])
                xparaminfo['start_angle'] = float(line[1])
                xparaminfo['osc_range'] = float(line[2])
                xparaminfo['gonio_coord'] = list(map(float, line[3:6]))
            elif ii == 2:
                line = allLines[ii].split()
                xparaminfo['wavelength'] = float(line[0])
                xparaminfo['incident_beam'] = list(map(float, line[1:4]))
            elif ii == 3:
                line = allLines[ii].split()
                xparaminfo['spacegroup'] = str(line[0])
                xparaminfo['unit_cell'] = list(map(float, line[1:7]))
            elif ii == 4:
                line = allLines[ii].split()
                xparaminfo['a_axis'] = list(map(float, line))
            elif ii == 5:
                line = allLines[ii].split()
                xparaminfo['b_axis'] = list(map(float, line))
            elif ii == 6:
                line = allLines[ii].split()
                xparaminfo['c_axis'] = list(map(float, line))
            elif ii == 7:
                line = allLines[ii].split()
                xparaminfo['NX'] = int(line[1])
                xparaminfo['NY'] = int(line[2])
                xparaminfo['pixelsize'] = float(line[3])
            elif ii == 8:
                line = allLines[ii].split()
                xparaminfo['beamX'] = float(line[0])
                xparaminfo['beamY'] = float(line[1])
                xparaminfo['detZ'] = float(line[2])
            elif ii == 9:
                line = allLines[ii].split()
                xparaminfo['det_axis1'] = list(map(float, line))
            elif ii == 10:
                xparaminfo['det_axis2'] = list(map(float, allLines[ii].split()))
            elif ii == 11:
                xparaminfo['det_axis3'] = list(map(float, allLines[ii].split()))
            else:
                pass
        return xparaminfo

    def collect_spots(self):
        spotfile = pathlib.Path(self.inData.get('XDSspotFile', ''))
        spot_list = []
        if not spotfile.exists():
            return spot_list
        with open(spotfile, 'r') as fh:
            for line in fh:
                x, y = map(lambda s: float(s), line.strip().split()[:2])
                spot_list.append((x, y))
        return spot_list

    def predict_miller_indices(self):
        crystal_symmetry = crystal.symmetry(unit_cell=self.xparamdict['unit_cell'],
                                            space_group_symbol=self.xparamdict['spacegroup'])
        d_min = self.inData.get('resolution', float('3'))
        d_max = None
        miller_set = miller.build_set(crystal_symmetry, anomalous_flag=True, d_min=d_min, d_max=d_max)
        return miller_set.indices()

    def calc_detector_projection(self):
        h_mat = self.predict_miller_indices()  # get predicted set of miller indices
        h_mat = np.array(h_mat)
        cell_mat = np.array([self.xparamdict['a_axis'], self.xparamdict['b_axis'], self.xparamdict['c_axis']]).transpose()
        rlp = np.linalg.inv(cell_mat)  # (a*, b*, c*) matrix

        # incident beam vector
        s0 = np.array(self.xparamdict['incident_beam'])
        beam_vector = s0 / np.linalg.norm(s0) / self.xparamdict['wavelength']

        rot_vector = np.array(self.xparamdict['gonio_coord'])  # obtained from rotation axis
        rot = np.empty((3, 3), dtype=np.float)
        rot[:, 1] = rot_vector / np.linalg.norm(rot_vector)
        rot[:, 0] = np.cross(rot[:, 1], beam_vector)
        rot[:, 0] /= np.linalg.norm(rot[:, 0])
        rot[:, 2] = np.cross(rot[:, 0], rot[:, 1])
        # detector coordinate system
        d_mat = np.array([self.xparamdict['det_axis1'], self.xparamdict['det_axis2'], self.xparamdict['det_axis3']])
        d_mat = d_mat.transpose()

        # 1st predicted spot as vector on Ewald sphere
        p0star = np.dot(h_mat, rlp)
        p0star_gonio = np.dot(p0star, rot)  # same predicted spot on gonio coordinate
        beam_gonio = np.dot(rot.transpose(), beam_vector)  # incident beam on gonio coordinate
        print(beam_gonio.shape)
        p0star_sq_dist = np.sum(p0star**2, axis=1)  # sq distance of predicted spot from the rotation axis

        pstar_gonio = np.empty(p0star_gonio.shape)  # a 3x3 matrix
        pstar_gonio[:, 2] = (-0.5*p0star_sq_dist - beam_gonio[1]*p0star_gonio[:, 1]) / beam_gonio[2]
        pstar_gonio[:, 1] = p0star_gonio[:, 1]
        pstar_gonio[:, 0] = p0star_sq_dist - pstar_gonio[:, 1]**2 - pstar_gonio[:, 2]**2  # sqrt taken at later stage

        sel_index = pstar_gonio[:, 0] > 0  # exclude blind region with this criteria, filters the matrices

        h_mat = h_mat[sel_index]
        p0star_gonio = p0star_gonio[sel_index]
        pstar_gonio = pstar_gonio[sel_index]

        pstar_gonio[:, 0] = np.sqrt(pstar_gonio[:, 0])

        projected_hkl = np.empty((0, 3), dtype=np.int)  # h,k,l as row vector
        position_on_detector = np.empty((0, 4), dtype=np.float)  # Find detector projection, x, y, phi, zeta

        for sign in (+1, -1):
            pstar_gonio[:, 0] *= sign
            phi = np.arctan2((p0star_gonio[:, 2]*pstar_gonio[:, 0] - p0star_gonio[:, 1]*pstar_gonio[:, 2]),
                             (p0star_gonio[:, 1]*pstar_gonio[:, 1] + p0star_gonio[:, 2]*pstar_gonio[:, 2]))

            refl_std = np.cross(pstar_gonio, beam_gonio)
            refl_std /= np.linalg.norm(refl_std, axis=1).reshape(refl_std.shape[0], 1)  # convert into column vector
            zeta = refl_std[:, 1]

            Svector = beam_vector + np.dot(pstar_gonio, rot.transpose())
            Svec_on_det = np.dot(Svector, d_mat)
            sel_index = self.xparamdict['detZ']*Svec_on_det[:, 2] > 0
            Svec_on_det, h_mat, phi, zeta = Svec_on_det[sel_index], h_mat[sel_index], phi[sel_index], zeta[sel_index]

            det_pos_x = self.xparamdict['beamX'] + self.xparamdict['detZ']*Svec_on_det[:, 0] / Svec_on_det[:, 2] \
                        / self.xparamdict['pixelsize']

            det_pos_y = self.xparamdict['beamY'] + self.xparamdict['detZ'] * Svec_on_det[:, 1] / Svec_on_det[:, 2] \
                        / self.xparamdict['pixelsize']

            projected_hkl = np.row_stack([projected_hkl, h_mat])
            position_on_detector = np.row_stack([position_on_detector,
                                                 np.column_stack([det_pos_x, det_pos_y, phi, zeta])])

        return projected_hkl, position_on_detector

    def get_spot_prediction(self, esd=3.0):

        mosaicity = self.inData.get('mosaicity', float('0.1'))
        frame = self.inData.get('frame', int(1))
        phi_from_expt = self.xparamdict['start_angle'] + \
                            self.xparamdict['osc_range']*(frame - self.xparamdict['start_angle'] + 0.5)
        phi_from_expt, mosaicity, osc_range = np.deg2rad([phi_from_expt, mosaicity, self.xparamdict['osc_range']])

        hkl_calc, projection = self.calc_detector_projection()
        phi_diff = np.fmod(projection[:, 2] - phi_from_expt, 2.0*np.pi)
        # put all phi differences within -pi to +pi range..

        phi_diff[phi_diff < - np.pi] += 2.0 * np.pi
        phi_diff[phi_diff > np.pi] -= 2.0 * np.pi
        # choose range of index to visualize..
        sel_index = np.abs(phi_diff) < (osc_range / 2. + esd * mosaicity / np.abs(projection[:, 3]))

        return hkl_calc[sel_index], projection[sel_index]


if __name__ == '__main__':
    inData = {
        "XDSspotFile": "/Users/shbasu/work/APICA/1x/SPOT.XDS",
        "XDSxparamFile": "/Users/shbasu/work/APICA/1x/XPARM.XDS"
    }

    pp = ExecSpotPrediction(inData)
    # xparam = pp.digestxparamFile(inData=inData)
    pmiller, pdata = pp.get_spot_prediction()
    spots = pp.collect_spots()
    distance_limit_in_pixel = 20

    out = open("tr1.adx", "w")
    nmatches = pp.nearest_neighbor(np.array(spots), pdata, distance_limit_in_pixel, out)
    '''
    with open("tst.adx", 'w') as ofh:
        for (h, k, l), (x, y, phi, zeta) in zip(pmiller, pdata):
            ofh.write('%d %d %d %d %d\n' % (x, y, h, k, l))
    '''

