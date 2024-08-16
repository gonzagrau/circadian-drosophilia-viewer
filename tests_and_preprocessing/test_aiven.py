from dotenv import load_dotenv
import pymysql
import os

def main():
  load_dotenv()
  db_host = os.getenv('DB_HOST')
  db_user = os.getenv('DB_USER')
  db_password = os.getenv('DB_PASS')

  timeout = 10
  connection = pymysql.connect(
    charset="utf8mb4",
    connect_timeout=timeout,
    cursorclass=pymysql.cursors.DictCursor,
    db="defaultdb",
    host=db_host,
    password=db_password,
    read_timeout=timeout,
    port=10184,
    user=db_user,
    write_timeout=timeout,
  )
    
  try:
    cursor = connection.cursor()
    cursor.execute("INSERT INTO mytest (id) VALUES (1), (2)")
    cursor.execute("SELECT * FROM mytest")
    print(cursor.fetchall())
  finally:
    connection.close()


if __name__ == '__main__':
  main()