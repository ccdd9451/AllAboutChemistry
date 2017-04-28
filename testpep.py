#!/usr/bin/env python
# encoding: utf-8

import os
import unittest as t
from pathlib import Path

import pepscanner as ps
from xilio import dump, write

ps.DEBUG = True
ps.BASEDIR = Path('/home/kh621/test')
ps.MD_NSTLIM = 5E4

def setuptest():
    testdir = ps.BASEDIR
    os.system('rm -r /home/kh621/test')
    os.mkdir(testdir)
    os.mkdir(testdir / 'Analysis')
    write(testdir / 'Analysis/summary_avg.EPTOT',
        "    5020.000           0.8468 ")

class ScannerTester(t.TestCase):
    testdir = ps.BASEDIR

    @t.skip('test method not decided')
    def test_singlerun(self):
        # ps.singlerun()
        pass

    def test_shelf(self):
        def shelftester():
            with ps.shelf_with_locker() as shelf:
                return shelf['queued'].pop()

        shelfpos = self.testdir / 'pep_queue'

        try:
            os.remove(shelfpos)
        except Exception:
            pass

        self.assertRaises(FileNotFoundError, shelftester)

        dump(shelfpos, {'queued': ['test']})
        self.assertEqual('test', shelftester())

    def test_setupfiles(self):
        workpos = self.testdir / 'GGG'
        pseqs = 'GLY GLY GLY'
        ps.setupfiles(workpos, pseqs, os.environ)
        self.assertTrue(os.path.isfile(workpos / 'ambsc'))
        self.assertTrue(os.path.isfile(workpos / 'tlsc'))
        self.assertTrue(os.path.isfile(workpos / 'inpcrd'))
        self.assertTrue(os.path.isfile(workpos / 'prmtop'))

    def test_aminoacid_converts(self):
        self.assertEqual('GLY GLY GLY GLY GLY',
                         ps.convert_short_to_long('GGGGG'))

    def test_parseptot(self):
        self.assertAlmostEqual(0.8468, ps.parseptot(self.testdir))


if __name__ == '__main__':
    setuptest()
    t.main()
