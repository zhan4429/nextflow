More formally, you can create functions that are defined as first class objects.
```
square = { it * it }
```

The curly brackets around the expression `it * it` tells the script interpreter to treat this expression as code. The `it` identifier is an implicit variable that represents the value that is passed to the function when it is invoked.

Once compiled the function object is assigned to the variable square as any other variable assignments shown previously. Now we can do something like this:

```
println square(9)
```
and get the value 81.

This is not very interesting until we find that we can pass the function square as an argument to other functions or methods. Some built-in functions take a function like this as an argument. One example is the `collect` method on lists:

```
[ 1, 2, 3, 4 ].collect(square)
```

This expression says: Create an array with the values 1, 2, 3 and 4, then call its collect method, passing in the closure we defined above. The collect method runs through each item in the array, calls the closure on the item, then puts the result in a new array, resulting in:
```
[ 1, 4, 9, 16 ]
```

