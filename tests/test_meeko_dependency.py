from sidock_ubuntu.cli import doctor_status


def test_doctor_requires_openbabel_and_never_checks_meeko():
    data = doctor_status()
    names = {item["name"] for item in data["dependencies"]}
    assert "Open Babel" in names
    assert not any("meeko" in name.lower() for name in names)
