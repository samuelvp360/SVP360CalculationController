# -*- coding: utf-8 -*-

import transaction
import os
from ZODB import FileStorage, DB
from BTrees.IOBTree import IOBTree


class MyZODB(object):

    def __init__(self):

        try:
            self.storage = FileStorage.FileStorage('DB/molecules_db.fs')
        except:
            os.remove('DB/molecules_db.fs.lock')
            self.storage = FileStorage.FileStorage('DB/molecules_db.fs')
        self.db = DB(self.storage)
        self.db.pack()
        self.connection = self.db.open()
        self.dbroot = self.connection.root()
        elements = ('config', 'molecules', 'projects')
        # elements_2 = ('jobs')
        # for element in elements_2:
        if not self.dbroot.get('jobs'):
            self.dbroot['jobs'] = IOBTree()
            self.commit()
        for element in elements:
            if not self.dbroot.get(element):
                self.dbroot[element] = {}
                self.commit()

    def close(self):
        self.connection.close()
        self.db.close()
        self.storage.close()

    def commit(self):
        transaction.commit()

    def set(self, kind, key, element):
        self.dbroot[kind][key] = element
        self.dbroot._p_changed = True
        self.commit()

    def get(self, kind, key):
        return self.dbroot[kind].get(key)

    def set_config(self, config_data):
        self.dbroot['config'] = config_data
        self.dbroot._p_changed = True
        self.commit()

    def update(self, kind, key, data):
        self.dbroot[kind][key].__dict__.update(data)
        if kind == 'molecules':
            self.dbroot._p_changed = True
        self.commit()

    @property
    def get_molecules_db(self):
        return self.dbroot['molecules'].values()

    @property
    def get_projects_db(self):
        return self.dbroot['projects'].values()

    @property
    def get_jobs_db(self):
        return self.dbroot['jobs'].values()

    @property
    def get_job_id(self):
        return len(self.dbroot['jobs'])

    @property
    def get_config(self):
        return self.dbroot.get('config')

    def check(self, kind, key):
        return key in self.dbroot[kind]

    def remove(self, kind, key):
        del self.dbroot[kind][key]
        self.dbroot._p_changed = True

