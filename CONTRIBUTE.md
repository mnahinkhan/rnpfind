## Development
Enable recommended Git Hooks as follows:
```
git config --local core.hooksPath .githooks/
```
The above will run the following to ensure code consistency every time you
commit:
 - [black](https://github.com/psf/black)
 - [isort](https://github.com/PyCQA/isort)

Also use [fit-commit](https://github.com/m1foley/fit-commit) to ensure
consistent commit message style.


### Steps for running the website
Use Python 3.8 or higher

Create a file called `db.env` (choose a db, username, and password):
```bash
# Set PostgreSQL credential env variables
POSTGRES_DB=******
POSTGRES_USER=*******
POSTGRES_PASSWORD=*****************
```

*Note:* In the following commands, docker might need `sudo` prepended depending
on your setup.

Pull the Docker images and run them
```bash
make prod # Might need sudo
```

Now you should be able to use `rnpfind` by navigating to localhost:80 (or the
public IP of your machine)

Once done, stop and remove by typing Ctrl+C or using `docker-compose down`
