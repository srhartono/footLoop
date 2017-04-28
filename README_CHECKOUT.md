It depends on what you want to do when you checkout that commit. If all you're doing is checking it out so you can build or test that revision, then there's nothing wrong with working with a detached head. Just remember to check out an actual branch before you make any commits (git checkout master, for example), so that you don't create commits that are not included in any branch.

If, however, you want to make more commits starting from that point, you should create a branch. If you make commits that are not referenced by a branch, they can easily get lost, and will eventually be cleaned up by git's garbage collector, as nothing refers to them. You can create a new branch by running:

git checkout -b newbranch ea3d5ed
To help visualize, here's are some diagrams demonstrating how working on a detached head differs from working on a branch.

Let's start out with 3 commits on master, A, B, and C. master is the current branch, so HEAD points to master, which points to commit C.

```
A  B  C
*--*--* <-- master <-- HEAD
```

Now if we commit, git will create a commit that has C as a parent (because that's the current commit, pointed to from HEAD via master), and will update master to point to that new commit. All of our commits are now in master, and HEAD points to the new commit through master.

```
A  B  C  D
*--*--*--* <-- master <-- HEAD
```

Now let's check out B, giving us a detached HEAD.

```
A  B  C  D
*--*--*--* <-- master
   ^
    \-- HEAD
```
Everything works fine here; we can look at all of the files, build our program, test it, etc. We can even create new commits; but if we do so, there's no branch that we're on, so we can't point any branch at that new commit. The only thing pointing at it is HEAD:

```
A  B  C  D
*--*--*--* <-- master
    \
     * <-- HEAD
     E
```

If we later decide to check out master again, there will be nothing referring to E.

```
A  B  C  D
*--*--*--* <-- master <-- HEAD
    \
     *
     E
```

Since there's nothing referring to it, it can be hard to find, and git considers commits with no references to be abandoned (they happen quite commonly if you rebase, or squash patches in, or do other fun history manipulation; they usually represent abandoned patches that you no longer care about). After a certain amount of time, git will consider it garbage, to be discarded the next time garbage collection runs.

So, instead of checking out a bare revision and getting a detached head, if you feel like you are going to make more commits, you should use git checkout -b branch B to create a branch and check it out. Now your commits won't be lost, as they will be included in a branch, that you can easily refer to, and merge later on.

```
A  B  C  D
*--*--*--* <-- master
   ^
    \-- branch <-- HEAD
```

If you forget to do this, and create commits off a branch, there's no need to worry. You can create a branch referring to the head revision with git checkout -b branch. If you have already switched back to the master branch, and realize that you forgot a stray commit, you can find it using git reflog, which will show you a history of what commits HEAD has pointed to over the last few days. Anything that's still in the reflog will not be garbage collected, and generally references are kept in the reflog for at least 30 days.

