## Nonmodifying Algorithms

for_each() - This is possibly the most important algorithm in this section, as it allows any unary function (i.e. a function of one argument) to be applied to each element in a range/container. Note that this function can actually also be modifying (hence why it is included below). It is often better to use a more specific algorithm, if one exists, than to use this, as specialist implementations will be more efficient.
count() - This returns the number of elements in a range or container.
count_if() - This counts how many elements in a range or container much a particular criterion.
min_element() - Returns the element that has the smallest value, making use of the < relation to perform comparison. It can accept a custom binary function to perform the comparison instead.
max_element() - Returns the element that has the largest value, making use of the > relation to perform comparison. It can accept a custom binary function to perform the comparison instead.
find() - Finds the first element in a range or container that equals a passed value.
find_if() - Finds the first element in a range or container that matches a particular criterion, rather than a passed value.
search_n() - This is like find AND find_if except that it looks for the first  occurances of such a value OR the first  occurances where a relational predicate is met.
search() - This searches for the first occurance of a subrange within a range/container and can do so either by first/last value of the subrange or via a predicate matching all the values of the desired first/last subrange.
find_end() - Similar to search, except that it finds the last occurance of such a subrange.
find_first_of() - Finds the first element in a subrange of a range or container. Can make use of a binary predicate function, otherwise uses direct value.
adjacent_find() - Returns the location (an iterator) to the first matching consecutive pair of values in a range/container. Can also match via a binary predicate.
equal() - Compares two ranges to see if they are equal.
mismatch() - Compares two ranges and returns a pair of iterators containing the points at which the ranges differ.
lexicographical_compare() - The lexicographical comparison is used to sort elements in a manner similar to how words are ordered in a dictionary. It can either use operator< or make use of a binary predicate function to perform the comparison.
Modifying Algorithms
for_each() - This is the same as the for_each we discussed above, but I have included it in the Modifying section to reinforce that it can be used this way too!
copy() - Copies a range/container of elements into another range.
copy_backward() - Copy a range/container of elements into another range, starting from the last element and working backwards.
transform() - Transform is quite a flexible algorithm. It works in two ways. A unary operation can be applied to the source range, on a per element basis, which ouputs the results in the destination range. A binary operation can be applied to both elements in the source and destination range, subsequently overwriting elements in the destination range.
merge() - Merge is intended to take two sorted ranges and combine them to produce a merged sorted range. However, it is possible to utilise unsorted ranges as arguments but this then leads to an unsorted merge! For this reason I've decided to include it in the modifying category, rather than the category for sorted ranges, below.
swap_ranges() - This swaps the elements of two ranges.
fill() - This replaces each element in a range with a specific value.
fill_n() - Similar to fill, but replaces the first  elements in a range with a specific value.
generate() - This replaces each element in a range with the result of an operation of a generator function.
generate_n() - Similar to generate, but replaces the first  elements in a range with the result of an operation of a generator function.
replace() - This replaces elements matching a specific value with another specific value.
replace_if() - This replaces elements matching a specific criterion with another specific value.
replace_copy() - Similar to replace, except that the result is copied into another range.
replace_copy_if() - Similar to replace_if, except that the result is copied into another range.

## Removal Algorithms

remove() - This removes elements from a range that match a specific value.
remove_if() - This removes elements from a range that match a specific criterion, as determined via a unary predicate.
remove_copy() - Similar to remove, except that elements are copied into another range.
remove_copy_if() - Similar to remove_if, except that elements are copied into another range.
unique() - This is quite a useful algorithm. It removes adjacent duplicate elements, i.e. consecutive elements with specific values.
unique_copy() - Similar to unique, except that it copies the elements into another range.
Mutating Algorithms
reverse() - This simply reverses the order of the elements in a range or container.
reverse_copy() - Similar to reverse, except that the results of the reversal are copied into another range.
rotate() - By choosing a 'middle' element in a range, this algorithm will cyclically rotate the elements such that the middle element becomes the first.
rotate_copy() - Similar to rotate, except that the result is copied into another range.
next_permutation() - This rearranges the elements in a range to produce the next lexicographically higher permutation, using operator<. It is also possible to use a binary predicate comparison function instead of operator<.
prev_permutation() - Similar to next_permutation, except that it rearranges to produce the next lexicographically lower permutation.
random_shuffle() - Rearranges the list of elements in a range in a random fashion. The source of randomness can be supplied as a random number generator argument.
partition() - Rearranges a range/container such that the elements matching a predicate are at the front. Does NOT guarantee relative ordering from the original range.
stable_partition() - Similar to partition, except that it does guarantee relative ordering from the original range.

## Sorting Algorithms
sort() - Sorts the elements into ascending order, using operator< or another supplied comparison function.
stable_sort() - This is similar to sort. It is used when you need to ensure that elements remain in the same order when they are "tied" for the same position. This often comes up when dealing with priorities of tasks. Note also that the performance guarantee is different.
partial_sort() - Similar to sort, except that it only sorts the first  elements and terminates after they're sorted.
partial_sort_copy() - Similar to partial_sort except that it copies the results into a new range.
nth_element() - This allows you to ensure that an element at position  is in the correct position, were the rest of the list to be sorted. It only guarantees that the elements preceeding  are less than in value (in the sense of operator<) and that proceeding elements are greater than in value.
partition() - See above for partition.
stable_partition() - See above for stable_partition.
make_heap() - Rearranges the elements in a range such that they form a heap, i.e. allowing fast retrieval of elements of the highest value and fast insertion of new elements.
push_heap() - This adds an element to a heap.
pop_heap() - This removes an element from a heap.
sort_heap() - This sorts the elements in a heap, with the caveat that the range is no longer a heap subsequent to the function call.

## Sorted Range Algorithms
binary_search() - Searches the range for any matches of a specific value.
includes() - Determines whether each element in one range is also an element in another range.
lower_bound() - Searches for the first element in a range which does not compare less than a specific value. Can also use a custom comparison function.
upper_bound() - Searches for the first element in a range which does not compare greater than a specific value. Can also use a custom comparison function.
equal_range() - Finds a subrange within a range which contains values equal to a specific value.
merge() - See above for merge.
set_union() - Creates a set union of two sorted ranges. Thus the destination range will include elements from either or both of the source ranges. Since this is a set operation, duplicates are eliminated.
set_intersection() - Creates a set intersection of two sorted ranges. Thus the destination range will include elements that exist only in both of the source ranges.
set_difference() - Creates a set difference of two sorted ranges. Thus the destation range will include elements from the first range that are not in the second range.
set_symmetric_difference() - This is the symmetric version of set_difference. It creates a set symmetric difference, which is formed by the elements that are found in one of the sets, but not in the other.
inplace_merge() - This combines two sorted ranges into a destination range, which is also sorted. It is also stable as it will preserve ordering for subranges.