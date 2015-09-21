function t = DB_gettable(DBconn, tablename)
% t = DB_gettable(DBconn, tablename)

t = fetch(DBconn, sprintf('select * from %s',tablename));
t = TableToCategorical(t,0);
