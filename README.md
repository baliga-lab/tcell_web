# T-Cell Mtb Infection Web Application

## Dependencies

 * Flask
 * Werkzeug
 * Jinja2
 * MySQLdb
 * numpy
 * pandas
 * rpy2

Install with

```
sudo pip install Flask
sudo pip install Werkzeug
sudo pip install Jinja2
sudo apt-get install python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose
sudo pip install rpy2
```

## Database setup
### Open mysql
mysql -u <user_name> -p
### Make the database
CREATE DATABASE tcell;
### Create structure of database
cat data/createDatabase.sql | mysql -u <user_name> -p
### Enter all data (sql file for database needs to be downloaded from website).
zcat ../tcell_web_data/tcell.sql.gz | mysql -u <user_name> -p

see settings.cfg, create MySQL database and adjust user/database settings
accordingly

## Run in Development mode

TCELL_SETTINGS=settings.cfg ./app.py
