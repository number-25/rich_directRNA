Command
`docker run number25/transigner:1.1.3 transigner --version`

```stderr
Traceback (most recent call last):
  File "/opt/conda/bin/transigner", line 5, in <module>
    from transigner.run_transigner import main
  File "/opt/conda/lib/python3.9/site-packages/transigner/run_transigner.py", line 4, in <module>
    from transigner import align, pre, em
  File "/opt/conda/lib/python3.9/site-packages/transigner/pre.py", line 158
    out_lns_pbase.append(f'{tname},{np.array2string(per_base_cov[i], separator=' ')}')
                                                                                     ^
SyntaxError: f-string: unmatched '('
```

Solution: Change the single quotes '' to double quotes ''

Error resolved.
