# Create documentation for the first time in a gh-pages branch
```shell
# install sphinx
pip install python-sphinx -y

# create docs folder
mkdir docs

# go in it and start sphinx
cd docs
sphinx-quickstart

# edit `conf.py` and `index.rst` to select options and create modules

# run sphinx API
sphinx-apidoc -o . ..

# make html page, this generates the _build directory
make html

# move built html directory outside the package to update gh-pages branch
mv _build/html ~/Downloads/

# remove the _build directory
rm -r _build

# go to repo root
cd ..

# checkout the gh-pages branch
git checkout gh-pages

# copy the contents of the built website
mv ~/Downloads/html/* .

# update to github
git add .
git commit -m "udpated documentation"
git push
```

# Update documentation in a gh-pages branch
```shell
# update files, `conf.py` and `index.rst`

# go to docs directory
cd docs

# run sphinx API
sphinx-apidoc -o . ..

# make html page
make html

# move built html directory outside the package to update gh-pages branch
mv _build/html ~/Downloads/

# remove the _build directory
rm -r _build

# go to repo root
cd ..

# checkout the gh-pages branch
git checkout gh-pages

# remove everythin non-hidden
rm -r *

# copy the contents of the built website
mv ~/Downloads/html/* .

# update to github
git add .
git commit -m "udpated documentation"
git push
```