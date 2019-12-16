group_by_evtno(table::TypedTables.Table) = TypedTables.Table(consgroupedview(table[sortperm(table.evtno)].evtno, Tables.columns(table)))
