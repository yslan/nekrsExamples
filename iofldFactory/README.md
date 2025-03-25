### A simple example of iofldFactory

v24pre (repo/next: Oct. 19 2025)

This cloned function demonstrates how to clone or costumize a iofldFactory so user can use it to dump arbitrary fields with arbitrary options.


See the original function `nrs_t::writeCheckpoint` (`nekrs/src/nrs/nrs.cpp`).

- Append user list of fields
  ```
  for (const auto &entry : userCheckpointFields) {
    checkpointWriter->addVariable(entry.first, entry.second);
  }
  ```

- Filter out unwanted elements
  ```
  auto elementFilter = [&]() 
  {
    std::vector<int> elementFilter;
    for(int e = 0; e < mesh->Nelements; e++) {
       auto zmaxLocal = std::numeric_limits<dfloat>::lowest();
       for(int i = 0; i < mesh->Np; i++) zmaxLocal = std::max(z[i + e * mesh->Np], zmaxLocal);
       if (zmaxLocal > zRecycLayer) elementFilter.push_back(e);
    }
    return elementFilter;
  }();
  iofld->writeElementFilter(elementFilter);
  ```
