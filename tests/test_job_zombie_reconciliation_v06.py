from sidock_ubuntu import job_store


def test_process_alive_rejects_zombie(monkeypatch):
    monkeypatch.setattr(job_store, 'process_identity', lambda pid: '12345')
    monkeypatch.setattr(job_store, '_linux_process_state', lambda pid: 'Z')
    assert job_store.process_alive(999, '12345') is False


def test_process_alive_accepts_matching_non_zombie(monkeypatch):
    monkeypatch.setattr(job_store, 'process_identity', lambda pid: '12345')
    monkeypatch.setattr(job_store, '_linux_process_state', lambda pid: 'S')
    assert job_store.process_alive(999, '12345') is True


def test_process_alive_rejects_pid_reuse(monkeypatch):
    monkeypatch.setattr(job_store, 'process_identity', lambda pid: 'new-start')
    monkeypatch.setattr(job_store, '_linux_process_state', lambda pid: 'R')
    assert job_store.process_alive(999, 'old-start') is False
