## Opening files
To access and work with files, use the `file` method, which returns a file system object given a file path string:
```
myFile = file('some/path/to/my_file.file')
```
**The file method can reference either files or directories, depending on what the string path refers to in the file system**.

When using the wildcard characters `*`, `?`, `[]` and `{}`, the argument is interpreted as a glob path matcher and the file method returns a list object holding the paths of files whose names match the specified pattern, or an empty list if no match is found:
```
listOfFiles = file('some/path/*.fa')
```

Two asterisks `(**)` in a glob pattern works like `*` but also searches through subdirectories.

By default, wildcard characters do not match directories or hidden files. For example, if you want to include hidden files in the result list, add the optional parameter hidden:

```
listWithHidden = file('some/path/*.fa', hidden: true)
```

