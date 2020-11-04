#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from ZODB import FileStorage, DB
import transaction


class MyZODB(object):

    def __init__(self):

        self.storage = FileStorage.FileStorage('DB/moleculesDB.fs')
        self.db = DB(self.storage)
        self.connection = self.db.open()
        self.dbroot = self.connection.root()

    def Close(self):

        self.connection.close()
        self.db.close()
        self.storage.close()

    def StoreMolecule(self, moleculeToStore):

        self.dbroot[moleculeToStore.GetInchikey] = moleculeToStore
        self.dbroot[moleculeToStore.GetInchikey].stored = True
        transaction.commit()

    def FetchMolecule(self, dbMolecule):

        if self.dbroot.get(dbMolecule) is not None:
            return self.dbroot[dbMolecule]
        else:
            return False

    def ShowAll(self):

        molObjectList = []
        molNames = []
        molKeys = []
        for key, molObject in self.dbroot.items():
            molObjectList.append(molObject)
            molNames.append(molObject.GetName)
            molKeys.append(molObject.GetInchikey)

        return molObjectList, molNames, molKeys

    def SearchMolecule(self, key):

        if self.dbroot.get(key) is not None:
            return True
        else:
            return False
