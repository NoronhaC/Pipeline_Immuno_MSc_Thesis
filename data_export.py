import mysql.connector


# Functions

def query_database(query):
    res = None
    cursor.execute(query)
    if query.lower().startswith('select') or query.lower().startswith('show'):
        res = cursor.fetchall()
    return res


# ---------------------------------------------------------------------------------
# Execution starts here

conn = mysql.connector.connect(
    host='localhost',
    user='root',
    password='root',
    database='biomedicaldb'
)

cursor = conn.cursor()

tables = [table[0] for table in query_database('SHOW TABLES')]
print(tables)

for table in tables:
    entries = query_database(f'SELECT * FROM biomedicaldb.{table}')
    print(f'Table Name: {table}')
    print(f'Entries: {entries}')
    print(f'Column Names: {cursor.column_names}')


query = "INSERT INTO fruits (name, quantity, country) VALUES (%s, %s, %s)"

val = [
  ("Banana", 40, "Mexico"),
  ("Mango", 15, "India"),
  ("Avocado", 37, "Mexico")
]

try:
    cursor.executemany(query, val)
    print(cursor.rowcount, "records inserted.")
except:
    print('Something went wrong.')

conn.commit()
