import sqlite3

conn=sqlite3.connect('data.sqlite')
c=conn.cursor()
c.execute("CREATE TABLE Friends(Id INTEGER PRIMARY KEY, Name TEXT);")
