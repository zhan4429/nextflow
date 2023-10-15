### pythonwork.nf
'''
nextflow.enable.dsl=2
process helloworld {

"""
#!/usr/bin/python

print("hello world")
x=1
y=5
print(x+y)
"""
}

workflow { 
    helloworld()
}
'''


```
nextflow run pythonwork.nf -process.echo
```
