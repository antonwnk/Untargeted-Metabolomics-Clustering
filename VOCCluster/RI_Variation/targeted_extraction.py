#!/usr/bin/env python2
# encoding: UTF-8
"""
    - This class will upload the dataset and targeted VOCs files which are located in the data folder.
    - It will search in the dataset about all of the targeted VOCs.
    - It uses threashold = 0.98 unless a different value has been passed.
    - At the end, a report will be generated and located in the data folder.
    @author: Yaser Alkhalifah, Jan - 2019
"""
import os
import numpy as np
import cosine_calculations


class Targeted:
    _sample_idx = 0
    _voc_num_idx = 1
    _ri_idx = 2
    _ri_end_idx = 3

    def __init__(self, data_file, targ_file, path, threshold=0.98):
        self.data_file = path + '/' + data_file
        self.targ_file = path + '/' + targ_file
        self.path = path
        self.threshold = threshold
        print 'Used threshold = ', self.threshold
        self.loaded_samples = []
        self.clusters = []
        self.retention_index = []
        self.ri_vary = []
        self.samples = []
        self.epsilon = 1
        # Load dataset, which includes heading and values
        f = open('{0}'.format(self.data_file), 'r')
        f = list(f)
        number_of_vocs = sum(1 for _ in f) - 1  # -1 heading
        number_of_columns = len(f[0].split(','))
        self.first_m_z = int(f[0].split(',')[3])
        self.last_m_z = int(f[0].split(',')[-1])

        # Create a matrix with a shape of (number_of_vocs X number_of_columns) filled with zeros.
        self.dataset = np.zeros((number_of_vocs,
                                 number_of_columns))
        # Fill in the matrix with the values.
        print 'Now uploading the dataset File.....'
        for line in range(1, len(f)):
            if int(float(f[line].strip().split(',')[self._sample_idx])) not in self.loaded_samples:
                self.loaded_samples.append(int(float(f[line].strip().split(',')[self._sample_idx])))
            for column in range(number_of_columns):
                self.dataset[line - 1][column] = int(float(f[line].strip().split(',')[column]))
        # Load Targeted VOCs into the targeted_vocs array
        print 'Number of samples that are loaded = ', len(self.loaded_samples)
        print 'dataset includes ', number_of_vocs, 'VOCs in all samples '
        print 'dataset includes ', number_of_columns, ' Columns, ', 'm/z values start from ', self.first_m_z, 'and end ', self.last_m_z
        print 'Now uploading the targeted VOCs File.....'
        f = open(self.targ_file, 'r')
        f = list(f)
        number_of_targeted = sum(1 for _ in f)
        print 'Number of targeted VOCs that have been uploaded : ', number_of_targeted
        self.targeted_vocs = np.zeros((number_of_targeted, 4))
        for line in range(number_of_targeted):
            for column in range(4):
                self.targeted_vocs[line][column] = int(float(f[line].strip().split(',')[column]))

    ##################################################################
    # This function to start clustering targeted compounds and stared in label(contains VOC index + cosine distance
    # value). Then calculate Delta RI for each segment. The segmented RI points and Delta RI are then stored in
    # self.ri_vary.
    def extract_targeted(self):
        """
            This function to extract all targeted VOCs from all samples.
        """
        print 'Searching about targeted VOCs ...'
        cosine = cosine_calculations.Cosine(self.dataset)
        # For each targeted VOC
        for target_idx, target_row in enumerate(self.targeted_vocs):

            # Find the targeted VOC
            for voc_idx, voc_row in enumerate(self.dataset):
                if target_row[self._sample_idx] == voc_row[self._sample_idx] and \
                        target_row[self._voc_num_idx] == voc_row[self._voc_num_idx]:
       b             self.clusters.append([voc_idx])
                    self.samples.append([target_row[self._sample_idx]])
                    self.retention_index.append([int(target_row[self._ri_idx])])
                    break

            # Extract all similar VOCs from all samples and cluster them...
            passes = 0
            while True:
                print 'Targeted VOC ', target_idx, ' pass number ', passes
                passes += 1
                no_update = True
                # For each VOC
                for voc_idx, voc_row in enumerate(self.dataset):
                    local_targeted_compound = []
                    # If found VOC that is in different sample and within RI variation range
                    if voc_row[self._sample_idx] not in self.samples[target_idx] and \
                            voc_row[self._ri_idx] in range(int(target_row[self._ri_idx]),
                                                           int(target_row[self._ri_end_idx])):
                        # For that VOC and all others after it
                        for local_index in range(voc_idx, len(self.dataset)):  # extract all VOCs in particular sample
                            # and cluster the highest cosine VOC
                            # If new VOC still within allowable RI range
                            if self.dataset[local_index][self._ri_idx] in range(int(target_row[self._ri_idx]),
                                                                                int(target_row[self._ri_end_idx])):
                                compared_m_z = cosine.normalisation(local_index)
                                # For each previously clustered VOC
                                for clu in self.clusters[target_idx]:
                                    targ_m_z = cosine.normalisation(clu)
                                    # If new and old VOCs similar enough store idx and distance of new
                                    if cosine.cosine_similarity(targ_m_z, compared_m_z) >= self.threshold:
                                        local_targeted_compound.append(
                                            (local_index, cosine.cosine_similarity(targ_m_z, compared_m_z)))
                            # If we go out of allowable range stop
                            else:
                                break
                    # If we found something in the allowable range of the currently explored sample
                    if len(local_targeted_compound) >= 1:
                        no_update = False
                        most_similar = max(local_targeted_compound, key=lambda tup: tup[1])[0]
                        self.clusters[target_idx].append(most_similar)
                        self.samples[target_idx].append(self.dataset[most_similar][self._sample_idx])
                        self.retention_index[target_idx].append(int(self.dataset[most_similar][self._ri_idx]))
                # When no more VOCs are added to the cluster
                if no_update:
                    n_samples = len(self.clusters[target_idx])
                    ri = self.retention_index[target_idx]
                    print 'VOC ', target_idx, ' was found in ', n_samples, ' samples.'
                    print 'Max RI = ', max(ri), ', min RI = ', min(ri), '.'
                    print 'Delta RI = ', max(ri) - min(ri)
                    break

        # Now testing all the clusters to find the max distance between 2 VOCs from the same cluster,and fill in self.ri_vary.
        for clu in range(len(self.clusters)):
            if clu == len(self.clusters) - 1:
                self.ri_vary.append([(max(self.retention_index[clu - 1]) + min(self.retention_index[clu])) / 2,
                                     max(self.retention_index[clu]) - min(self.retention_index[clu])])
            else:
                self.ri_vary.append([(max(self.retention_index[clu]) + min(self.retention_index[clu + 1])) / 2,
                                     max(self.retention_index[clu]) - min(self.retention_index[clu])])

            for in_clust in range(len(self.clusters[clu])):
                for in_clust2 in range(len(self.clusters[clu])):
                    if in_clust != in_clust2:
                        if cosine.cosine_similarity(cosine.normalisation(self.clusters[clu][in_clust]),
                                                    cosine.normalisation(self.clusters[clu][in_clust2])) < self.epsilon:
                            self.epsilon = cosine.cosine_similarity(cosine.normalisation(self.clusters[clu][in_clust]),
                                                                    cosine.normalisation(self.clusters[clu][in_clust2]))
        print 'Epsilon = ', self.epsilon
        self.ri_vary.sort(key=lambda tup: tup[0], reverse=False)
        self.ri_variation_report()

    # Create a CSV report contains all Targeted clusters...
    def ri_variation_report(self):
        """
            All targeted VOCs that have been found will be reported in
            a CSV file in order to check the alignment.
            """
        result_file = open(self.path + '/' + 'RI_Report.csv', 'w')
        result_file.write('Sample#' + ",")
        result_file.write('VOC#' + ",")
        result_file.write('RI' + ",")
        for m_z in range(self.first_m_z, self.last_m_z + 1):
            result_file.write(str(m_z) + ",")
        result_file.write("\n")
        for clu in range(len(self.clusters)):
            for in_clust in range(len(self.clusters[clu])):
                result_file.write(str(self.dataset[int(self.clusters[clu][in_clust])][0]) + ",")
                result_file.write(str(self.dataset[int(self.clusters[clu][in_clust])][1]) + ",")
                result_file.write(str(self.dataset[int(self.clusters[clu][in_clust])][2]) + ",")
                for voc in range(3, len(self.dataset[self.clusters[clu][in_clust]])):
                    result_file.write(str(self.dataset[int(self.clusters[clu][in_clust])][voc]) + ",")
                result_file.write("\n")
            result_file.write('##########################################################')
            result_file.write("\n")
            result_file.write('' + ",")
            result_file.write('MaxRI' + ",")
            result_file.write(str(max(self.retention_index[clu])) + ",")
            result_file.write("\n")
            result_file.write('' + ",")
            result_file.write('MinRI' + ",")
            result_file.write(str(min(self.retention_index[clu])) + ",")
            result_file.write("\n")
            result_file.write('' + ",")
            result_file.write('Delta RI' + ",")
            result_file.write(str(max(self.retention_index[clu]) - min(self.retention_index[clu])) + ",")
            result_file.write("\n")
            result_file.write('' + ",")
            result_file.write('RI COV %' + ",")
            result_file.write(str(np.std(self.retention_index[clu]) / np.mean(self.retention_index[clu]) * 100) + ",")
            result_file.write("\n")
            result_file.write('##########################################################')
            result_file.write("\n")
        result_file.close()
        print 'Report file has been generated and saved in the data folder.'
