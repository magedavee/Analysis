import sqlite3 as lite
import sys

def setupTables(name='P50Events.sqlite'):
    con=lite.connect(name)
    with con:
            
        cur = con.cursor()    
        cur.execute("CREATE TABLE P50Muon(event INT,pulseTop REAL,pulseBot REAL,renormTop REAL,renormBot REAL,length REAL,time REAL)")
def insertEvent(event,pulseTop,pulseBot,renormTop,renormBot,length,time,name='P50Events.sqlite',table='P50Muon') :
    con=lite.connect(name)
    with con:
        cmd="INSERT INTO %s VALUES(%d,%f,%f,%f,%f,%f,%f)" %(table,event,pulseTop,pulseBot,renormTop,renormBot,length,time)
        # print cmd
        cur = con.cursor()    
        cur.execute(cmd)
# def addCol(name,Type,table='P50Muon'):
        # cmd="ALTER TABLE %s ADD COLUMN %s %s " %(table,name,Type)
        # print cmd
