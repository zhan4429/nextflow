### hello.nf
```
nextflow.enable.dsl=2
process sayHello {

"""
echo 'Hello world!' > /cluster/home/yzhang85/notes/nextflow/hello_world.txt
"""
}

workflow {

    sayHello()
}
```

```
$ nextflow run hello.nf
```


### hello2.nf
```
nextflow.enable.dsl=2
process sayHello {

"""
echo 'Hello world!'
"""
}

workflow {

    sayHello()
}
```

```
$ nextflow run hello2.nf
N E X T F L O W  ~  version 21.10.6
Launching `hello2.nf` [voluminous_tuckerman] - revision: 9a35148624
executor >  local (1)
[7f/b0b4c4] process > sayHello [100%] 1 of 1 ✔
```

```
$ nextflow run hello2.nf -process.echo
N E X T F L O W  ~  version 21.10.6
Launching `hello2.nf` [extravagant_saha] - revision: 9a35148624
executor >  local (1)
[94/4b436e] process > sayHello [100%] 1 of 1 ✔
Hello world!
```

**Note:**  If we do not add `-process.hello`, the echo information will not be shown.

### hellopython.nf
```
nextflow.enable.dsl=2
process sayHello {

"""
#!/usr/bin/python3
print("Hello world")
"""

}

workflow{
    sayHello()

}
```

```
$ nextflow run hellopython.nf -process.echo
N E X T F L O W  ~  version 21.10.6
Launching `hellopython.nf` [ecstatic_archimedes] - revision: 619aa7b0e5
executor >  local (1)
[0d/e46566] process > sayHello [100%] 1 of 1 ✔
Hello world
```
