# -*- coding: utf-8 -*-

import transaction
import os
from ZODB import FileStorage, DB


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
        elements = ('config', 'molecules')
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

    def set(self, molecule):
        self.dbroot['molecules'][molecule.inchi_key] = molecule
        self.dbroot._p_changed = True
        self.commit()

    def set_config(self, config_data):
        self.dbroot['config'] = config_data
        self.dbroot._p_changed = True
        self.commit()

    def update(self, key, data):
        self.dbroot['molecules'][key].__dict__.update(data)
        self.dbroot._p_changed = True
        self.commit()

    @property
    def get_db(self):
        return self.dbroot['molecules'].values()

    def get_molecule(self, key):
        return self.dbroot['molecules'].get(key)

    @property
    def get_config(self):
        return self.dbroot.get('config')

    def check(self, key):
        return key in self.dbroot['molecules']

    def remove(self, key):
        del self.dbroot['molecules'][key]
        self.dbroot._p_changed = True
        self.commit()

