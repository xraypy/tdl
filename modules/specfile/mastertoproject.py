'''
Master to Project Parser
Author: Craig Biwer (cbiwer@uchicago.edu)
4/20/2012
'''

import h5py
import sys

INTEGRATION_PARAMETERS = {'bgrflag': 1,
                          'cnbgr': 5,
                          'colormap': 'None',
                          'compress': 1,
                          'cpow': 2,
                          'ctan': False,
                          'cwidth': 15,
                          'filter': False,
                          'geom': 'psic',
                          'integrated': False,
                          'nline': 1,
                          'rnbgr': 5,
                          'roi': [],
                          'rpow': 0,
                          'rtan': False,
                          'rwidth': 15}
CORRECTION_PARAMETERS = {'bad_pixel_map': 'None',
                         'bad_point': False,
                         'image_changed': True,
                         'image_max': -1,
                         'pixel_map_changed': True,
                         'real_image_max': -1,
                         'rotangle': 0,
                         'sample_angles': [],
                         'sample_diameter': 10,
                         'sample_polygon': [],
                         'scale': 1.e6}
DETECTOR_PARAMETERS = {'beam_slits': [],
                       'det_slits': [],
                       'name': 'pilatus'
                      }
                      
def master_to_project(master_file, desired_scans, project_file):
    '''Convert a list of scans from a master file
        (potentially with some listed attributes)
        to a project file.
        
        Inputs:
            -master_file: the full directory path to the master file
            -desired_scans: a nested dictionary describing what scans
                            to pull from the master file and write to
                            the project file. Format is:
                                {spec_file: {scan: {attribute: value}}}
            -project_file: the full directory path of the desired
                            output file
                            
    '''
    
    read_this = h5py.File(master_file, 'r')
    write_this = h5py.File(project_file, 'a')
    
    # The starting point number (this only works if
    # the file consists of sequentially numbered
    # points, starting at 1 with no gaps).
    point_counter = len(write_this) + 1
    
    for spec_name in desired_scans.keys():
        for scan_number in desired_scans[spec_name].keys():
            read_head = read_this[spec_name][scan_number]
            num_points = len(read_head['point_data'])
            #add_time = 0
            for i in range(num_points):
                point_name = '%6.6i' % point_counter
                point_group = write_this.create_group(point_name)
                point_counter += 1
                
                specified_attrs = desired_scans[spec_name][scan_number]
                
                # The epoch of the point (specfile creation time + point offset)
                file_epoch = read_head.attrs.get('init_epoch', 0)
                epoch_loc = list(read_head['point_labs']).index('Epoch')
                point_epoch = read_head['point_data'][i][epoch_loc]
                
                # The unique identifier
                point_group.attrs['name'] = str(spec_name + ':S' + \
                                                scan_number +':P' + str(i) + \
                                                '/' + str(num_points) + ':' + \
                                                str(point_epoch + file_epoch))
                # The scan type
                point_group.attrs['type'] = read_head.attrs.get('s_type', 'NA')
                # Info about the point; initialized to the spec command
                point_group.attrs['info'] = read_head.attrs.get('cmd', 'NA')
                # The geometry of the point
                point_group.attrs['geom'] = \
                        specified_attrs.get('geom', INTEGRATION_PARAMETERS['geom'])
                # A comment field; initialized to 'Initialized'
                point_group.attrs['hist.1'] = 'Initialized'
                # The time the point was taken (saved as an epoch)
                # NOTE: This is determined by adding the Epoch offset of the
                # point to the epoch at the beginning of the file
                point_group.attrs['date_stamp'] = point_epoch + file_epoch
                #The energy
                point_group.attrs['energy'] = read_head.attrs.get('energy', 'NA')
                
                # The angles
                ang_labels = ['chi', 'del', 'eta', 'mu', 'nu', 'phi']
                point_group.create_dataset('angle_labels', data=ang_labels)
                pang_labels = ['chi', 'TwoTheta', 'theta', 'Psi', 'Nu', 'phi']
                ang_values = []
                for j in range(6):
                    ang_lbl = ang_labels[j]
                    pang_lbl = pang_labels[j]
                    if ang_lbl in read_head['point_labs']:
                        ang_pos = list(read_head['point_labs']).index(ang_lbl)
                        ang_val = read_head['point_data'][i][ang_pos]
                    elif pang_lbl in read_head['point_labs']:
                        ang_pos = list(read_head['point_labs']).index(pang_lbl)
                        ang_val = read_head['point_data'][i][ang_pos]
                    else:
                        ang_pos = list(read_head['param_labs']).index(pang_lbl)
                        ang_val = read_head['param_data'][ang_pos]
                    ang_values.append(ang_val)
                point_group.create_dataset('angle_values', data=ang_values)
                
                # The real-space lattice followed by the recip-space lattice
                lattice_start = list(read_head['param_labs']).index('g_aa')
                #lattice_stop = list(read_head['param_labs']).index('g_ga_s')
                ltc_lbls = ['real_a', 'real_b', 'real_c', 'real_alpha',
                            'real_beta', 'real_gamma', 'recip_a', 'recip_b',
                            'recip_c', 'recip_alpha', 'recip_beta',
                            'recip_gamma']
                point_group.create_dataset('lattice_labels', data=ltc_lbls)
                point_group.create_dataset('lattice_values',
                                           data=read_head['param_data']\
                                           [lattice_start:lattice_start+12])
                
                # Q
                h_loc = list(read_head['point_labs']).index('H')
                k_loc = list(read_head['point_labs']).index('K')
                L_loc = list(read_head['point_labs']).index('L')
                h_val = read_head['point_data'][i][h_loc]
                k_val = read_head['point_data'][i][k_loc]
                L_val = read_head['point_data'][i][L_loc]
                point_group.create_dataset('Q', data=[h_val, k_val, L_val])
                
                # The or's (all of or0 followed by all of or1)
                # HKLs
                or_start = list(read_head['param_labs']).index('g_h0')
                or_zero = list(read_head['param_data'][or_start:or_start+3])
                or_one = list(read_head['param_data'][or_start+3:or_start+6])
                # Angles
                or_start += 6
                or_zero.extend(read_head['param_data'][or_start:or_start+6])
                or_one.extend(read_head['param_data'][or_start+6:or_start+12])
                # Lambdas
                or_start += 12
                or_zero.extend(read_head['param_data'][or_start:or_start+1])
                or_one.extend(read_head['param_data'][or_start+1:or_start+2])
                or_zero.extend(or_one)
                or_labs = ['or0_h', 'or0_k', 'or0_L', 'or0_del', 'or0_eta',
                           'or0_chi', 'or0_phi', 'or0_nu', 'or0_mu',
                           'or0_lambda', 'or1_h', 'or1_k', 'or1_L', 'or1_del',
                           'or1_eta', 'or1_chi', 'or1_phi', 'or1_nu', 'or1_mu',
                           'or1_lambda']
                point_group.create_dataset('or_labels', data=or_labs)
                point_group.create_dataset('or_values', data=or_zero)
                
                # Azimuth vector
                haz_start = list(read_head['param_labs']).index('g_haz')
                haz_data = read_head['param_data'][haz_start:haz_start+3]
                point_group.create_dataset('haz', data=haz_data)
                
                # The position values
                split_index = list(read_head['point_labs']).index('Epoch')
                pos_lbls = read_head['point_labs'][:split_index]
                point_group.create_dataset('position_labels', data=pos_lbls)
                pos_values = read_head['point_data'][i][:split_index]
                point_group.create_dataset('position_values', data=pos_values)
                
                # The scaler values
                sclr_lbls = read_head['point_labs'][split_index:]
                point_group.create_dataset('scaler_labels', data=sclr_lbls)
                sclr_values = read_head['point_data'][i][split_index:]
                point_group.create_dataset('scaler_values', data=sclr_values)
                
                # The detector
                det_group = point_group.create_group('det_0')
                # Name
                no_show = DETECTOR_PARAMETERS.get('name', 'NA')
                det_group.attrs['name'] = specified_attrs.get('name', no_show)
                # Data
                det_group.attrs['data'] = read_head['point_data'][i]
                # Image, if it exists
                try:
                    det_group.create_dataset('image_data',
                                             data=read_head['image_data'][i],
                                             compression='szip')
                except:
                    pass
                # Integration parameters
                int_labels = ['bgrflag', 'cnbrg', 'compress', 'cpow', 'ctan',
                              'cwidth', 'filter', 'integrated', 'nline',
                              'rnbgr', 'roi', 'rpow', 'rtan', 'rwidth']
                det_group.create_dataset('int_labels', data=int_labels)
                int_values = []
                for label in int_labels:
                    no_show = INTEGRATION_PARAMETERS.get(label, 'NA')
                    int_values.append(str(specified_attrs.get(label, no_show)))
                det_group.create_dataset('int_values', data=int_values)
                # Correction parameters
                corr_labels = ['bad_pixel_map', 'bad_point',
                               'image_changed', 'image_max',
                               'pixel_map_changed', 'real_image_max',
                               'rotangle', 'sample_angles', 'sample_diameter',
                               'sample_polygon', 'scale']
                det_group.create_dataset('corr_labels', data=corr_labels)
                corr_values = []
                for label in corr_labels:
                    no_show = CORRECTION_PARAMETERS.get(label, 'NA')
                    corr_values.append(str(specified_attrs.get(label, no_show)))
                det_group.create_dataset('corr_values', data=corr_values)
                # Detector parameters
                det_labels = ['beam_slits', 'det_slits']
                det_group.create_dataset('det_labels', data=det_labels)
                det_values = []
                for label in det_labels:
                    no_show = DETECTOR_PARAMETERS.get(label, 'NA')
                    det_values.append(str(specified_attrs.get(label, no_show)))
                det_group.create_dataset('det_values', data=det_values)
                # Results
                res_labels = ['alpha', 'beta', 'ctot', 'F', 'F_changed', 'Ferr',
                              'I', 'I_c', 'I_r', 'Ibgr', 'Ibgr_c', 'Ibgr_r',
                              'Ierr', 'Ierr_c', 'Ierr_r']
                res_values = [0, 0, 0, 0, True, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                det_group.create_dataset('result_labels', data=res_labels)
                det_group.create_dataset('result_values', data=res_values)
                
                '''
                for attr in DEFAULT_ATTRIBUTES.keys():
                    if attr not in desired_scans[spec_name][scan_number].keys():
                        point_group.attrs[attr] = DEFAULT_ATTRIBUTES[attr]
                    else:
                        point_group.attrs[attr] = str(desired_scans[spec_name][scan_number][attr])
                point_group.create_dataset('raw_data', data=read_head['point_data'][i])
                try:
                    point_group.create_dataset('image_data', data=read_head['image_data'][i], compression='szip')
                except:
                    pass'''
    
    read_this.close()
    write_this.close()