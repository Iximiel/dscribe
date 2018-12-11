# Build docs, copy to correct docs folder, delete build
cd docs/src
sphinx-apidoc -o ./doc ../../dscribe
make html
cp build/html/ ../
rm -r build


# Push changes to docs
git config --global user.email "travis@travis-ci.org"
git config --global user.name "Travis CI"
git add ./docs
git commit -m "Travis documentation build: $TRAVIS_BUILD_NUMBER"
git push --quiet https://SINGROUP:$GH_TOKEN@github.com/SINGROUP/dscribe development &>/dev/null
