function t = DB_tableheaders(DBconn, tablename)
%t = DB_tableheaders(DBconn, tablename)

t = fetch(DBconn, sprintf('pragma table_info(''%s'')',tablename));
