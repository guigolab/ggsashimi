# Development

## Make a release

Check the current development release:
```
git checkout dev
$ ./ggsashimi.py --version
ggsashimi v0.7.0-dev
```

### Release the current development version

Create a 0.7.0 release from the development version above:

```
git checkout dev
git checkout -b release/0.7.0
```

Use [`bump2version`](https://github.com/c4urself/bump2version) to make the release:

```
bump2version release
```

This command strips the `-dev` suffix from the release and create a `v0.7.0` tag.

### Major release

Create a 1.0.0 release from the development version above:

```
git checkout dev
git checkout -b release/1.0.0
```

Use [`bump2version`](https://github.com/c4urself/bump2version) to update the version and create the tag:

```
bum2version --no-tag major
bump2version release
```

The first command bumps the version to `1.0.0-dev` and skips the creation of a tag. The second command strips the `-dev` suffix and create the tag.

