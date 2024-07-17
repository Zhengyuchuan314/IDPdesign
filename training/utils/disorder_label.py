import time
import random
import os
import multiprocessing
from multiprocessing import JoinableQueue
from iupred2a.iupred2a_lib import iupred
from Bio import SeqIO
import torch
import json
import numpy as np


def get_directories(id):
    real_id = id.split('UniRef50_')[1]

    # 一级目录分为三个，A0A，UPI，short，二级目录采用id的第2-4位
    if real_id.startswith('A0A'):
        primary_dir = 'A0A'
        secondary_dir = real_id[4:7]
    elif real_id.startswith('UPI'):
        primary_dir = 'UPI'
        secondary_dir = real_id[4:7]
    else:
        primary_dir = 'short'
        secondary_dir = real_id[1:4]

    return primary_dir, secondary_dir, real_id


class GetSeq(multiprocessing.Process):

    def __init__(self, t_name, queue, file_path):
        multiprocessing.Process.__init__(self, name=t_name)
        self.queue = queue  # 任务队列
        self.t_name = t_name
        self.file_path = file_path

    def run(self):
        for seq in SeqIO.parse(self.file_path, "fasta"):
            seq_data = {'ID': seq.id, 'seq': seq.seq}
            self.queue.put(seq_data)
        print('all sequences have been loaded')
        self.queue.join()
        print("all sequences have been get")


class DisorderCalculator(multiprocessing.Process):

    def __init__(self, t_name, queue, output_dir):
        multiprocessing.Process.__init__(self, name=t_name)
        self.queue = queue		    # task queue
        self.t_name = t_name
        self.output_dir = output_dir

    def run(self):
        while True:
            seq_data = self.queue.get()
            if seq_data is None:
                break
            else:
                try:
                    disorder = iupred(seq_data['seq'],'short')[0]
                    primary_dir, secondary_dir, real_id = get_directories(seq_data['ID'])
                    if not os.path.exists(os.path.join(self.output_dir, primary_dir)):
                        os.mkdir(os.path.join(self.output_dir, primary_dir))
                    if not os.path.exists(os.path.join(self.output_dir, primary_dir, secondary_dir)):
                        os.mkdir(os.path.join(self.output_dir, primary_dir, secondary_dir))
                    data = {'id': seq_data['ID'], 'seq':str(seq_data['seq']), 'disorder': disorder}
                    #如果要保存为pt格式
                    # torch.save(data, os.path.join(self.output_dir, primary_dir, secondary_dir, f'{real_id}.pt'))
                    with open(os.path.join(self.output_dir, primary_dir, secondary_dir, f'{real_id}.json'),'w') as f:
                        json.dump(data, f)
                    self.queue.task_done()
                except:
                    self.queue.task_done()


if __name__ == '__main__':
    workers = 12

    q = JoinableQueue(maxsize=24)

    # seq_getter
    producer  = GetSeq("get_seq", q, 'uniref50_first_50000_lines.fasta')
    producer.start()

    # calculator
    consumer_list = []
    for i in range(workers):
        consumer = DisorderCalculator(f"calculator {i}", q, './')
        # consumer.daemon = True
        consumer.start()
        consumer_list.append(consumer)

    producer.join()

    for i in range(workers):
        q.put(None)
    for consumer in consumer_list:
        consumer.join()


